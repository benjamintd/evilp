#include "dive.h"
#include "common.h"
#include "exact.h"
#include "native_tableau.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(c) do { if (!(c)) { fprintf(stderr,"dive driver line %d: %s\n",__LINE__,#c); return 1; } } while (0)
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

int main(int argc, char **argv) {
    if (argc != 4) return 64;
    const char *mode=argv[3];
    LabModel model={0};CHECK(!lab_model_read(&model,argv[1]));
    LabLpOptions lp={.time_limit=30,.iteration_limit=100000,.max_dimension=16000000};
    LabLpSession *main_session=lab_lp_session_create(&model,&lp);CHECK(main_session);
    LabLpResult root={0};lab_lp_session_solve(main_session,&root);
    if(strcmp(root.status,"optimal")||!root.has_primal){
        printf("{\"root_status\":\"%s\",\"status\":\"not_attempted\"}\n",root.status);
        lab_lp_result_free(&root);lab_lp_session_free(main_session);lab_model_free(&model);return 0;
    }
    size_t n=model.cols;double *lower=lab_alloc(n,sizeof(double)),*upper=lab_alloc(n,sizeof(double));
    memcpy(lower,model.col_lower,n*sizeof(double));memcpy(upper,model.col_upper,n*sizeof(double));
    LabDiveOptions options=lab_dive_options(1);
    if(!strcmp(mode,"cold"))options.cold=1;
    if(!strcmp(mode,"no-backtracks"))options.max_backtracks=0;
    if(!strcmp(mode,"one-step"))options.max_steps=1;
    if(!strcmp(mode,"tiny-budget"))options.time_limit=1e-100;
    if(!strcmp(mode,"row-cap"))options.max_rows=0;
    if(!strcmp(mode,"invalid-root")){
        CHECK(n>0&&isfinite(lower[0]));lower[0]-=1;
    }
    if(!strcmp(mode,"lp-abstain"))lp.iteration_limit=1;
    double *saved_lower=lab_alloc(n,sizeof(double)),*saved_upper=lab_alloc(n,sizeof(double));
    double *saved_point=lab_alloc(n,sizeof(double));
    memcpy(saved_lower,lower,n*sizeof(double));memcpy(saved_upper,upper,n*sizeof(double));
    memcpy(saved_point,root.primal,n*sizeof(double));
    uint64_t original=fingerprint(&model);
    LabNativeLimits limits={.max_rows=4,.max_terms=1024,.max_basis_dimension=1024,
        .max_basis_nnz=100000,.max_canonical_nnz=100000,.max_bits=8192,
        .max_total_bits=16000000,.deadline=lab_now()+1};
    LabNativeBatch before={0},after={0};
    LabNativeStatus native_before=lab_lp_session_export_native(main_session,&limits,&before);
    LabDiveResult result={0};lab_primal_dive(&model,&model,lower,upper,root.primal,&lp,&options,&result);
    CHECK(original==fingerprint(&model));CHECK(!memcmp(saved_point,root.primal,n*sizeof(double)));
    CHECK(!memcmp(saved_lower,lower,n*sizeof(double)));CHECK(!memcmp(saved_upper,upper,n*sizeof(double)));
    limits.deadline=lab_now()+1;
    LabNativeStatus native_after=lab_lp_session_export_native(main_session,&limits,&after);
    CHECK(native_before==native_after);
    if(native_before==LAB_NATIVE_OK){
        CHECK(before.dimension==after.dimension&&before.count==after.count);
        CHECK(!memcmp(before.basis,after.basis,before.dimension*sizeof(size_t)));
        for(size_t i=0;i<before.dimension;i++)CHECK(!mpq_cmp(before.basic_values[i],after.basic_values[i]));
    }
    LabLpResult repeated={0};lab_lp_session_solve(main_session,&repeated);
    CHECK(!strcmp(repeated.status,"optimal")&&repeated.has_primal);
    CHECK(!memcmp(saved_point,repeated.primal,n*sizeof(double)));
    if(result.primal){
        FILE*f=fopen(argv[2],"wx");CHECK(f);fprintf(f,"=obj= %.17g\n",result.objective);
        for(size_t j=0;j<n;j++)fprintf(f,"%s %.17g\n",model.col_name[j],result.primal[j]);
        CHECK(!ferror(f)&&!fclose(f));
    }
    fputs("{\"root_status\":\"optimal\",\"status\":",stdout);lab_json_string(stdout,lab_dive_status_name(result.status));
    fputs(",\"message\":",stdout);lab_json_string(stdout,result.message);
    printf(",\"candidate\":%s,\"objective\":%.17g,\"steps\":%zu,\"backtracks\":%zu,"
           "\"lp_solves\":%zu,\"lp_iterations\":%zu,\"warm_starts\":%zu,\"rounded_checks\":%zu,"
           "\"final_resolves\":%zu,\"elapsed\":%.17g,\"inputs_unchanged\":true,"
           "\"main_session_unchanged\":true,\"native_before\":%d,\"native_after\":%d}\n",
           result.primal?"true":"false",result.objective,result.steps,result.backtracks,result.lp_solves,
           result.lp_iterations,result.warm_starts,result.rounded_checks,result.final_resolves,result.elapsed,
           (int)native_before,(int)native_after);
    lab_dive_result_free(&result);CHECK(!result.primal);lab_dive_result_free(&result);
    lab_native_batch_free(&before);lab_native_batch_free(&after);
    lab_lp_result_free(&root);lab_lp_result_free(&repeated);lab_lp_session_free(main_session);
    free(lower);free(upper);free(saved_lower);free(saved_upper);free(saved_point);lab_model_free(&model);
    return 0;
}
