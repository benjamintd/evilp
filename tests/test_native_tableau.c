#include "native_tableau.h"
#include "common.h"
#include "exact.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) do { if (!(c)) { fprintf(stderr,"tableau contract line %d: %s\n",__LINE__,#c); return -1; } } while (0)
static size_t exports, fault_rejections, stale_rejections;
static FILE *evidence;
static LabLpOptions options(void) {
    return (LabLpOptions){.time_limit=120,.iteration_limit=10000,.max_dimension=4096};
}
static LabNativeLimits limits(void) {
    return (LabNativeLimits){.max_rows=4,.max_terms=1024,.max_basis_dimension=1024,
        .max_basis_nnz=100000,.max_canonical_nnz=100000,.max_bits=8192,
        .max_total_bits=16000000,.deadline=lab_now()+20};
}
static char *name(const char *prefix,size_t i) {
    char text[64]; snprintf(text,sizeof text,"%s%zu",prefix,i); return lab_strdup(text);
}
static void make_model(LabModel *m,size_t rows,size_t cols,const double *a,
                       const double *rl,const double *ru,const double *lo,const double *hi,
                       const double *obj,const unsigned char *integer) {
    memset(m,0,sizeof *m);m->rows=rows;m->cols=cols;
    m->row_name=lab_alloc(rows,sizeof(char*));m->col_name=lab_alloc(cols,sizeof(char*));
    m->row_lower=lab_alloc(rows,sizeof(double));m->row_upper=lab_alloc(rows,sizeof(double));
    m->col_lower=lab_alloc(cols,sizeof(double));m->col_upper=lab_alloc(cols,sizeof(double));
    m->objective=lab_alloc(cols,sizeof(double));m->integer=lab_alloc(cols,sizeof(unsigned char));
    m->col_start=lab_alloc(cols+1,sizeof(size_t));m->row_index=lab_alloc(rows*cols,sizeof(size_t));
    m->value=lab_alloc(rows*cols,sizeof(double));
    for(size_t i=0;i<rows;i++){m->row_name[i]=name("r",i);m->row_lower[i]=rl[i];m->row_upper[i]=ru[i];}
    for(size_t j=0;j<cols;j++){
        m->col_name[j]=name("x",j);m->col_lower[j]=lo[j];m->col_upper[j]=hi[j];
        m->objective[j]=obj[j];m->integer[j]=integer[j];m->col_start[j]=m->nnz;
        for(size_t i=0;i<rows;i++)if(a[i*cols+j]!=0){m->row_index[m->nnz]=i;m->value[m->nnz++]=a[i*cols+j];}
    }
    m->col_start[cols]=m->nnz;
    m->envelope=lab_source_envelope_new(rows,cols,m->nnz);
    for(size_t k=0;k<m->nnz;k++)m->envelope->value[k]=(LabInterval){m->value[k],m->value[k]};
    for(size_t i=0;i<rows;i++){
        m->envelope->row_lower[i]=(LabInterval){rl[i],rl[i]};m->envelope->row_upper[i]=(LabInterval){ru[i],ru[i]};
    }
    for(size_t j=0;j<cols;j++){
        m->envelope->col_lower[j]=(LabInterval){lo[j],lo[j]};m->envelope->col_upper[j]=(LabInterval){hi[j],hi[j]};
        m->envelope->col_lower_source[j]=lo[j];m->envelope->col_upper_source[j]=hi[j];
    }
    m->envelope->complete=1;
}
static void promote_exact(LabModel *m) {
    m->exact=lab_exact_model_new(m->rows,m->cols,m->nnz);m->parsed_decimal=1;
    for(size_t k=0;k<m->nnz;k++)mpq_set_d(m->exact->value[k],m->value[k]);
    for(size_t i=0;i<m->rows;i++){
        if(isfinite(m->row_lower[i]))mpq_set_d(m->exact->row_lower[i],m->row_lower[i]);
        if(isfinite(m->row_upper[i]))mpq_set_d(m->exact->row_upper[i],m->row_upper[i]);
    }
    for(size_t j=0;j<m->cols;j++){
        if(isfinite(m->col_lower[j]))mpq_set_d(m->exact->col_lower[j],m->col_lower[j]);
        if(isfinite(m->col_upper[j]))mpq_set_d(m->exact->col_upper[j],m->col_upper[j]);
        mpq_set_d(m->exact->objective[j],m->objective[j]);
        m->exact->col_lower_source[j]=m->col_lower[j];m->exact->col_upper_source[j]=m->col_upper[j];
    }
}
static uint64_t bytes(uint64_t h,const void *memory,size_t n) {
    const unsigned char *p=memory;for(size_t i=0;i<n;i++){h^=p[i];h*=UINT64_C(1099511628211);}return h;
}
static uint64_t qhash(uint64_t h,const mpq_t q){char *s=mpq_get_str(NULL,10,q);h=bytes(h,s,strlen(s)+1);free(s);return h;}
static uint64_t fingerprint(const LabModel *m) {
    uint64_t h=UINT64_C(1469598103934665603);
    h=bytes(h,m->value,m->nnz*sizeof(double));h=bytes(h,m->row_index,m->nnz*sizeof(size_t));
    h=bytes(h,m->col_start,(m->cols+1)*sizeof(size_t));h=bytes(h,m->row_lower,m->rows*sizeof(double));
    h=bytes(h,m->row_upper,m->rows*sizeof(double));h=bytes(h,m->col_lower,m->cols*sizeof(double));
    h=bytes(h,m->col_upper,m->cols*sizeof(double));h=bytes(h,m->objective,m->cols*sizeof(double));
    h=bytes(h,m->integer,m->cols*sizeof(unsigned char));
    for(size_t i=0;i<m->rows;i++)h=bytes(h,m->row_name[i],strlen(m->row_name[i])+1);
    for(size_t j=0;j<m->cols;j++)h=bytes(h,m->col_name[j],strlen(m->col_name[j])+1);
    if(m->exact){
        for(size_t k=0;k<m->nnz;k++)h=qhash(h,m->exact->value[k]);
        for(size_t i=0;i<m->rows;i++){h=qhash(h,m->exact->row_lower[i]);h=qhash(h,m->exact->row_upper[i]);}
        for(size_t j=0;j<m->cols;j++){h=qhash(h,m->exact->col_lower[j]);h=qhash(h,m->exact->col_upper[j]);h=qhash(h,m->exact->objective[j]);}
        h=qhash(h,m->exact->offset);h=bytes(h,m->exact->col_lower_source,m->cols*sizeof(double));h=bytes(h,m->exact->col_upper_source,m->cols*sizeof(double));
    }
    if(m->envelope){
        h=bytes(h,m->envelope->value,m->nnz*sizeof(LabInterval));
        h=bytes(h,m->envelope->row_lower,m->rows*sizeof(LabInterval));h=bytes(h,m->envelope->row_upper,m->rows*sizeof(LabInterval));
        h=bytes(h,m->envelope->col_lower,m->cols*sizeof(LabInterval));h=bytes(h,m->envelope->col_upper,m->cols*sizeof(LabInterval));
    }
    return h;
}
static void printq(FILE *f,const mpq_t q){char*s=mpq_get_str(NULL,10,q);fprintf(f,"\"%s\"",s);free(s);}
static void printmatrix(FILE *f,const LabNativeMatrix *m){
    fprintf(f,"{\"rows\":%zu,\"cols\":%zu,\"start\":[",m->rows,m->cols);
    for(size_t j=0;j<=m->cols;j++)fprintf(f,"%s%zu",j?",":"",m->start[j]);
    fputs("],\"row\":[",f);for(size_t k=0;k<m->nnz;k++)fprintf(f,"%s%zu",k?",":"",m->row[k]);
    fputs("],\"value\":[",f);for(size_t k=0;k<m->nnz;k++){if(k)fputc(',',f);printq(f,m->value[k]);}fputs("]}",f);
}
static void dump(const LabNativeBatch *b){
    if(!evidence)return;
    fprintf(evidence,"{\"scope\":\"ORIGINAL_MIP_IN_BOX; exporter algebra only\",\"dimension\":%zu,\"total_rows\":%zu,\"source\":",b->dimension,b->total_rows);
    printmatrix(evidence,&b->source_matrix);fputs(",\"basis\":",evidence);printmatrix(evidence,&b->basis_matrix);
    fputs(",\"canonical\":",evidence);printmatrix(evidence,&b->canonical_matrix);
    fputs(",\"rhs\":[",evidence);for(size_t i=0;i<b->dimension;i++){if(i)fputc(',',evidence);printq(evidence,b->rhs[i]);}
    fputs("],\"basic_values\":[",evidence);for(size_t i=0;i<b->dimension;i++){if(i)fputc(',',evidence);printq(evidence,b->basic_values[i]);}
    fputs("],\"basis_columns\":[",evidence);for(size_t i=0;i<b->dimension;i++)fprintf(evidence,"%s%zu",i?",":"",b->basis[i]);
    fputs("],\"integer\":[",evidence);for(size_t j=0;j<b->columns;j++)fprintf(evidence,"%s%d",j?",":"",(int)b->integer[j]);
    fputs("],\"row_id\":[",evidence);for(size_t i=0;i<b->dimension;i++)fprintf(evidence,"%s%zu",i?",":"",b->row_id[i]);
    fputs("],\"box\":[",evidence);for(size_t j=0;j<b->columns;j++){
        if(j)fputc(',',evidence);fputs("{\"lower\":",evidence);if(b->box[j].has_lower)printq(evidence,b->box[j].lower);else fputs("null",evidence);
        fputs(",\"upper\":",evidence);if(b->box[j].has_upper)printq(evidence,b->box[j].upper);else fputs("null",evidence);fputc('}',evidence);
    }
    fputs("],\"source_row_box\":[",evidence);for(size_t i=0;i<b->source_rows;i++){
        if(i)fputc(',',evidence);fputs("{\"lower\":",evidence);if(b->source_row_box[i].has_lower)printq(evidence,b->source_row_box[i].lower);else fputs("null",evidence);
        fputs(",\"upper\":",evidence);if(b->source_row_box[i].has_upper)printq(evidence,b->source_row_box[i].upper);else fputs("null",evidence);fputc('}',evidence);
    }
    fputs("],\"tableaux\":[",evidence);
    for(size_t r=0;r<b->count;r++){
        const LabNativeTableau*t=&b->tableau[r];if(r)fputc(',',evidence);
        fprintf(evidence,"{\"basic\":%zu,\"basis_position\":%zu,\"beta\":",t->basic_column,t->basis_position);printq(evidence,t->beta);
        fputs(",\"signed_basic_scale\":",evidence);printq(evidence,t->signed_basic_scale);
        fputs(",\"basic_shift\":",evidence);printq(evidence,t->basic_shift);
        fputs(",\"weights\":[",evidence);for(size_t i=0;i<t->weight_count;i++){if(i)fputc(',',evidence);printq(evidence,t->weights[i]);}
        fputs("],\"terms\":[",evidence);for(size_t k=0;k<t->terms;k++){
            const LabNativeTerm*q=&t->term[k];fprintf(evidence,"%s{\"kind\":%d,\"index\":%zu,\"upper\":%d,\"integer\":%d,\"canonical_column\":%zu,\"coefficient\":",k?",":"",(int)q->kind,q->index,q->upper,q->is_integer,q->canonical_column);printq(evidence,q->coefficient);
            fputs(",\"anchor\":",evidence);printq(evidence,q->anchor);fputs(",\"canonical_scale\":",evidence);printq(evidence,q->canonical_scale);fputs(",\"canonical_offset\":",evidence);printq(evidence,q->canonical_offset);fputc('}',evidence);
        }
        fputs("],\"equalities\":[",evidence);for(size_t k=0;k<t->equalities;k++){fprintf(evidence,"%s{\"index\":%zu,\"multiplier\":",k?",":"",t->equality[k].index);printq(evidence,t->equality[k].multiplier);fputc('}',evidence);}
        fputs("],\"fixed\":[",evidence);for(size_t k=0;k<t->fixed;k++){fprintf(evidence,"%s{\"index\":%zu,\"multiplier\":",k?",":"",t->fix[k].index);printq(evidence,t->fix[k].multiplier);fputc('}',evidence);}fputs("]}",evidence);
    }
    fputs("]}\n",evidence);fflush(evidence);
}
static int valid_batch(const LabNativeBatch *b){
    CHECK(b->count>0&&b->scope&&strstr(b->scope,"no independent"));
    for(size_t i=0;i<b->count;i++){
        const LabNativeTableau*t=&b->tableau[i];CHECK(t->basis_identity_verified&&t->native_identity_verified);
        CHECK(b->integer[t->basic_column]);
        for(size_t k=0;k<t->terms;k++){
            const LabNativeTerm*q=&t->term[k];CHECK(q->kind!=LAB_NATIVE_VAR||q->index!=t->basic_column);
            if(q->is_integer)CHECK(q->kind==LAB_NATIVE_VAR&&b->integer[q->index]&&mpz_cmp_ui(mpq_denref(q->anchor),1)==0);
        }
    }
    return 0;
}
static int export_checked(LabLpSession*s,const LabModel*m,const LabLpResult*r,LabNativeBatch*b){
    uint64_t h=fingerprint(m);LabLpResult saved=*r;LabNativeLimits l=limits();
    LabNativeStatus status=lab_lp_session_export_native(s,&l,b);
    CHECK(status==LAB_NATIVE_OK);CHECK(h==fingerprint(m)&&!memcmp(&saved,r,sizeof saved));
    CHECK(valid_batch(b)==0);dump(b);exports++;return 0;
}
static int scalar_profiles(void){
    const unsigned char integer[]={1};
    const double lower[]={0,2,-INFINITY,0,0},upper[]={16,18,-2,16,16};
    const double rl[]={-INFINITY,-INFINITY,-7,.5,.5},ru[]={3,7,INFINITY,3,3};
    const double cost[]={-1,-1,1,-1,1};
    const long numerator[]={3,7,-7,3,1},denominator[]={2,2,2,2,4};
    for(size_t k=0;k<5;k++)for(int dense=0;dense<2;dense++)for(int scale=0;scale<2;scale++){
        const double a[]={2};LabModel m;make_model(&m,1,1,a,&rl[k],&ru[k],&lower[k],&upper[k],&cost[k],integer);
        LabLpOptions o=options();o.dense=dense;o.no_scaling=!scale;LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);
        LabLpResult r;lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal")&&r.has_primal);
        LabNativeBatch b;lab_native_batch_init(&b);CHECK(export_checked(s,&m,&r,&b)==0);CHECK(b.count==1);
        CHECK(mpq_cmp_si(b.tableau[0].beta,numerator[k],(unsigned long)denominator[k])==0);
        if(scale&&k<2)CHECK(mpq_cmp_ui(b.tableau[0].signed_basic_scale,16,1)==0);
        if(k==2)CHECK(mpq_sgn(b.tableau[0].signed_basic_scale)<0);
        int row_term=0;for(size_t i=0;i<b.tableau[0].terms;i++)if(b.tableau[0].term[i].kind==LAB_NATIVE_ROW){row_term++;CHECK(b.tableau[0].term[i].upper==(k<2||k==3));}
        CHECK(row_term>0);lab_native_batch_free(&b);lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    }
    return 0;
}
static int shift_classification(void){
    const double a[]={2,2},rl[]={-INFINITY},lo0[]={0,.25},hi0[]={16,16.25},obj[]={-1,0};const unsigned char integer[]={1,1};
    for(int integral=0;integral<2;integral++){
        double lo[]={lo0[0],integral?2:lo0[1]},hi[]={hi0[0],integral?18:hi0[1]},ru[]={integral?7:3};
        LabModel m;make_model(&m,1,2,a,rl,ru,lo,hi,obj,integer);LabLpOptions o=options();LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);
        LabLpResult r;lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));LabNativeBatch b={0};CHECK(export_checked(s,&m,&r,&b)==0);
        int seen=0;for(size_t k=0;k<b.tableau[0].terms;k++){LabNativeTerm*t=&b.tableau[0].term[k];if(t->kind==LAB_NATIVE_VAR&&t->index==1){seen++;CHECK(t->is_integer==integral);}}
        CHECK(seen==1);lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    }
    return 0;
}
static int equality_fixed_and_typing(void){
    const unsigned char ix[]={1};const double lo[]={0},hi[]={16},obj[]={-1};
    const double a[]={2,4},side[]={1,2};LabModel m;make_model(&m,2,1,a,side,side,lo,hi,obj,ix);
    LabLpOptions o=options();LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);LabLpResult r;lab_lp_session_solve(s,&r);
    CHECK(!strcmp(r.status,"optimal")&&r.has_primal&&r.primal[0]==.5);
    LabNativeBatch b={0};CHECK(export_checked(s,&m,&r,&b)==0);
    CHECK(b.dimension<b.total_rows&&b.tableau[0].equalities>0);
    CHECK(b.tableau[0].terms==0&&mpq_cmp_ui(b.tableau[0].beta,1,2)==0);
    /* LP feasibility is established above. A later integer-only builder may
     * derive 0 >= 1/2; that must never be typed as original LP infeasibility. */
    lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    const double af[]={2,3},sf[]={1.75},lf[]={0,.25},uf[]={16,.25},of[]={-1,0};const unsigned char inf[]={1,0};
    make_model(&m,1,2,af,sf,sf,lf,uf,of,inf);promote_exact(&m);s=lab_lp_session_create(&m,&o);CHECK(s);lab_lp_session_solve(s,&r);
    CHECK(!strcmp(r.status,"optimal"));CHECK(export_checked(s,&m,&r,&b)==0);
    CHECK(b.tableau[0].equalities==1&&b.tableau[0].fixed==1&&b.tableau[0].fix[0].index==1);
    lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);return 0;
}
static int exact_decimal_rhs(void){
    const double a[]={1},rl[]={-INFINITY},ru[]={1.5},lo[]={0},hi[]={16},obj[]={-1};const unsigned char ix[]={1};
    LabModel m;make_model(&m,1,1,a,rl,ru,lo,hi,obj,ix);promote_exact(&m);
    CHECK(lab_exact_decimal(m.exact->value[0],"1.000000000000000001")==0);
    LabLpOptions o=options();LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);LabLpResult r;lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));
    LabNativeBatch b={0};CHECK(export_checked(s,&m,&r,&b)==0);mpq_t expected;mpq_init(expected);mpq_set_ui(expected,3,2);mpq_div(expected,expected,m.exact->value[0]);
    CHECK(mpq_cmp(expected,b.tableau[0].beta)==0&&mpq_cmp_ui(expected,3,2)!=0);
    CHECK(mpq_cmp(m.exact->value[0],b.source_matrix.value[0])==0);mpq_clear(expected);
    lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    CHECK(valid_batch(&b)==0);lab_native_batch_free(&b);return 0;
}
static int bound_slacks_and_deleted_row(void){
    const double a[]={2,2},rl[]={-INFINITY},ru[]={3},lo[]={0,0},hi[]={1,1},obj[]={-1,-1};const unsigned char ix[]={1,1};
    LabModel m;make_model(&m,1,2,a,rl,ru,lo,hi,obj,ix);LabLpOptions o=options();LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);
    LabLpResult r;lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));LabNativeBatch b={0};CHECK(export_checked(s,&m,&r,&b)==0);
    int bound_seen=0;for(size_t k=0;k<b.tableau[0].terms;k++)if(b.tableau[0].term[k].kind==LAB_NATIVE_BOUND){bound_seen++;CHECK(!b.tableau[0].term[k].is_integer);}
    CHECK(bound_seen>0);lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    const double aa[]={0,2},rlo[]={0,-INFINITY},rhi[]={0,3},l[]={0},u[]={16},c[]={-1};const unsigned char ii[]={1};
    make_model(&m,2,1,aa,rlo,rhi,l,u,c,ii);s=lab_lp_session_create(&m,&o);CHECK(s);lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));CHECK(export_checked(s,&m,&r,&b)==0);
    CHECK(b.dimension<b.total_rows&&b.row_id[0]!=0);int row_seen=0;
    for(size_t k=0;k<b.tableau[0].terms;k++)if(b.tableau[0].term[k].kind==LAB_NATIVE_ROW){row_seen++;CHECK(b.tableau[0].term[k].index==1);}
    CHECK(row_seen>0);lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);return 0;
}
static int fixed_decimal_shift(void){
    const double a[]={1,1},side[]={.6},lo[]={0,.1},hi[]={16,.1},obj[]={-1,0};const unsigned char ix[]={1,0};
    LabModel m;make_model(&m,1,2,a,side,side,lo,hi,obj,ix);promote_exact(&m);
    CHECK(!lab_exact_decimal(m.exact->row_lower[0],".6")&&!lab_exact_decimal(m.exact->row_upper[0],".6"));
    CHECK(!lab_exact_decimal(m.exact->col_lower[1],".1")&&!lab_exact_decimal(m.exact->col_upper[1],".1"));
    LabLpOptions o=options();LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);LabLpResult r;lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));
    LabNativeBatch b={0};CHECK(export_checked(s,&m,&r,&b)==0);CHECK(mpq_cmp_ui(b.tableau[0].beta,1,2)==0);
    CHECK(b.tableau[0].fixed==1&&mpq_cmp_ui(b.box[1].lower,1,10)==0);
    lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);return 0;
}
static int solve_discard_result(LabLpSession *s){LabLpResult r;lab_lp_session_solve(s,&r);int okay=!strcmp(r.status,"optimal");lab_lp_result_free(&r);return okay?0:-1;}
static int generations_and_failures(void){
    const double a[]={2},rl[]={-INFINITY},ru[]={3},lo[]={0},hi[]={16},obj[]={-1};const unsigned char ix[]={1};
    LabModel m;make_model(&m,1,1,a,rl,ru,lo,hi,obj,ix);LabLpOptions o=options();LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);
    LabNativeLimits l=limits();LabNativeBatch b={0};CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_STALE);stale_rejections++;
    LabLpResult r;lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));CHECK(export_checked(s,&m,&r,&b)==0);
    CHECK(lab_lp_session_set_bounds(s,0,0,16)==0);CHECK(export_checked(s,&m,&r,&b)==0);
    CHECK(lab_lp_session_set_bounds(s,0,-1,16)==-1);CHECK(export_checked(s,&m,&r,&b)==0);
    CHECK(lab_lp_session_set_bounds(s,0,1,16)==0);CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_STALE);stale_rejections++;
    lab_lp_result_free(&r);lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));CHECK(export_checked(s,&m,&r,&b)==0);
    CHECK(lab_lp_session_set_bounds(s,0,2,1)==0);CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_STALE);stale_rejections++;
    lab_lp_result_free(&r);lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"infeasible"));CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_STALE);stale_rejections++;
    CHECK(lab_lp_session_set_bounds(s,0,0,16)==0);lab_lp_result_free(&r);lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));CHECK(export_checked(s,&m,&r,&b)==0);
    LabNativeBatch old=b;uint64_t model_before=fingerprint(&m);LabLpResult saved=r;
    LabNativeLimits bad=l;bad.deadline=lab_now()-1;CHECK(lab_lp_session_export_native(s,&bad,&b)==LAB_NATIVE_DEADLINE);CHECK(!memcmp(&old,&b,sizeof b));
    bad=l;bad.max_terms=0;CHECK(lab_lp_session_export_native(s,&bad,&b)==LAB_NATIVE_RESOURCE);CHECK(!memcmp(&old,&b,sizeof b));
    bad=l;bad.max_bits=1;CHECK(lab_lp_session_export_native(s,&bad,&b)==LAB_NATIVE_RESOURCE);CHECK(!memcmp(&old,&b,sizeof b));
    bad=l;bad.max_total_bits=1;CHECK(lab_lp_session_export_native(s,&bad,&b)==LAB_NATIVE_RESOURCE);CHECK(!memcmp(&old,&b,sizeof b));
    bad=l;bad.max_basis_dimension=1;CHECK(lab_lp_session_export_native(s,&bad,&b)==LAB_NATIVE_RESOURCE);CHECK(!memcmp(&old,&b,sizeof b));
    bad=l;bad.max_basis_nnz=1;CHECK(lab_lp_session_export_native(s,&bad,&b)==LAB_NATIVE_RESOURCE);CHECK(!memcmp(&old,&b,sizeof b));
    bad=l;bad.max_canonical_nnz=1;CHECK(lab_lp_session_export_native(s,&bad,&b)==LAB_NATIVE_RESOURCE);CHECK(!memcmp(&old,&b,sizeof b));
    bad=l;bad.max_rows=0;CHECK(lab_lp_session_export_native(s,&bad,&b)==LAB_NATIVE_INVALID);CHECK(!memcmp(&old,&b,sizeof b));
    int succeeded=0;
    for(size_t k=0;k<100;k++){
        old=b;lab_native_test_fail_after(k);LabNativeStatus result=lab_lp_session_export_native(s,&l,&b);lab_native_test_fail_after(SIZE_MAX);
        if(result==LAB_NATIVE_OK){succeeded=1;break;}
        CHECK(result==LAB_NATIVE_RESOURCE);CHECK(!memcmp(&old,&b,sizeof b));CHECK(valid_batch(&b)==0);fault_rejections++;
    }
    CHECK(succeeded&&fault_rejections>20);CHECK(model_before==fingerprint(&m)&&!memcmp(&saved,&r,sizeof r));
    /* Work.result now refers to an expired stack object. Export must never
     * use it as scratch space or update its counters. */
    CHECK(solve_discard_result(s)==0);CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_OK);CHECK(valid_batch(&b)==0);
    lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);return 0;
}
static int abstentions_and_row_cap(void){
    const double a[]={2},rl[]={-INFINITY},ru[]={3},lo[]={0},hi[]={16},obj[]={-1};unsigned char integer[]={0};
    LabModel m;make_model(&m,1,1,a,rl,ru,lo,hi,obj,integer);LabLpOptions o=options();LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);LabLpResult r;lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));
    LabNativeLimits l=limits();LabNativeBatch b={0};CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_NO_ROWS);CHECK(!b.count&&!b.box);
    lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    integer[0]=1;make_model(&m,1,1,a,rl,ru,lo,hi,obj,integer);m.parsed_decimal=1;s=lab_lp_session_create(&m,&o);CHECK(s);lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_OK);lab_native_batch_free(&b);
    lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    /* A proved point enclosure is provenance; a parsed model with neither
     * exact rationals nor an enclosure must still abstain. */
    make_model(&m,1,1,a,rl,ru,lo,hi,obj,integer);m.parsed_decimal=1;
    lab_source_envelope_free(m.envelope);m.envelope=NULL;
    s=lab_lp_session_create(&m,&o);CHECK(s);lab_lp_session_solve(s,&r);
    CHECK(!strcmp(r.status,"optimal"));
    CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_SOURCE_UNAVAILABLE);
    lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    const double side[]={1},lf[]={-INFINITY},uf[]={INFINITY};make_model(&m,1,1,a,side,side,lf,uf,obj,integer);s=lab_lp_session_create(&m,&o);CHECK(s);lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_FREE_SPLIT);
    lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);
    const double aa[]={2,0,0,2},rrl[]={-INFINITY,-INFINITY},rru[]={3,3},ll[]={0,0},uu[]={16,16},oo[]={-1,-1};const unsigned char ii[]={1,1};
    make_model(&m,2,2,aa,rrl,rru,ll,uu,oo,ii);s=lab_lp_session_create(&m,&o);CHECK(s);lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));l.max_rows=1;CHECK(lab_lp_session_export_native(s,&l,&b)==LAB_NATIVE_OK);CHECK(b.count==1&&b.row_limit_reached&&b.tableau[0].basic_column==0);CHECK(valid_batch(&b)==0);dump(&b);exports++;
    lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);lab_model_free(&m);return 0;
}
int main(int argc,char **argv){
    if(argc==2){evidence=fopen(argv[1],"wx");if(!evidence)return 2;}
    if(scalar_profiles()||shift_classification()||equality_fixed_and_typing()||exact_decimal_rhs()||bound_slacks_and_deleted_row()||fixed_decimal_shift()||generations_and_failures()||abstentions_and_row_cap())return 1;
    if(evidence)fclose(evidence);
    printf("{\"passed\":true,\"checked_exports\":%zu,\"allocation_rejections\":%zu,\"stale_rejections\":%zu,\"cut_insertion\":false,\"independent_cut_validity_claim\":false}\n",exports,fault_rejections,stale_rejections);
    return 0;
}
