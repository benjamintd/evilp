#include "model.h"
#include "exact.h"
#include "root_cut.h"
#include "common.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"getter control %d: %s\n",__LINE__,#c);return 1;}}while(0)
int main(int argc,char **argv){
 CHECK(argc==3);LabModel m={0};CHECK(!lab_model_read(&m,argv[1]));CHECK(!m.exact&&m.parsed_decimal&&m.rows==200&&m.cols==300);
 /* These controls intentionally exercise point-only provenance. Native-only
  * rational coverage is validated separately against the original source. */
 lab_exact_model_free(m.native_exact);m.native_exact=NULL;
 mpq_t q;mpq_init(q);CHECK(!lab_source_exact_coefficient(&m,0,q)&&mpq_cmp_ui(q,2,1)==0);
 mpq_set_ui(q,42,1);m.envelope->complete=0;CHECK(lab_source_exact_coefficient(&m,0,q)&&mpq_cmp_ui(q,42,1)==0);m.envelope->complete=1;
 LabInterval saved=m.envelope->value[0];m.envelope->value[0].upper=nextafter(2,INFINITY);CHECK(lab_source_exact_coefficient(&m,0,q)&&mpq_cmp_ui(q,42,1)==0);m.envelope->value[0]=saved;
 LabExactModel *bad=lab_exact_model_new(m.rows+1,m.cols,m.nnz);m.exact=bad;CHECK(lab_source_exact_coefficient(&m,0,q)&&mpq_cmp_ui(q,42,1)==0);m.exact=NULL;lab_exact_model_free(bad);
 saved=m.envelope->col_lower[0];m.envelope->col_lower[0]=(LabInterval){nextafter(0,-INFINITY),nextafter(0,INFINITY)};
 CHECK(lab_source_exact_column_bound(&m,0,0,q)&&mpq_cmp_ui(q,42,1)==0);m.col_lower[0]=.25;
 CHECK(!lab_source_exact_column_bound(&m,0,0,q)&&mpq_cmp_ui(q,1,4)==0);mpq_set_ui(q,42,1);
 CHECK(lab_source_exact_column_value(&m,0,0,0,q)&&mpq_cmp_ui(q,42,1)==0);m.col_lower[0]=0;m.envelope->col_lower[0]=saved;
 saved=m.envelope->row_upper[199];m.envelope->row_upper[199].lower=nextafter(m.row_upper[199],-INFINITY);
 CHECK(lab_source_exact_row_bound(&m,199,1,q)&&mpq_cmp_ui(q,42,1)==0);
 double *lo=malloc(m.rows*sizeof(double)),*hi=malloc(m.rows*sizeof(double));CHECK(lo&&hi);memcpy(lo,m.row_lower,m.rows*sizeof(double));memcpy(hi,m.row_upper,m.rows*sizeof(double));lo[199]=-INFINITY;hi[199]=INFINITY;LabModel view=m;view.row_lower=lo;view.row_upper=hi;
 LabLpOptions o={.time_limit=30,.iteration_limit=10000,.max_dimension=4096};LabLpSession *session=lab_lp_session_create(&view,&o);CHECK(session);LabLpResult r;lab_lp_session_solve(session,&r);CHECK(!strcmp(r.status,"optimal"));
 LabRootCutLimits limits=lab_root_cut_limits(lab_now()+10);LabRootCutBatch b;lab_root_cut_batch_init(&b);CHECK(lab_root_cut_derive(session,&m,&limits,&b)==LAB_ROOT_CUT_OK&&b.count==1);
 FILE *f=fopen(argv[2],"wx");CHECK(f&&fputs(b.cut[0].certificate,f)>=0&&!fclose(f));
 lab_root_cut_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(session);free(lo);free(hi);m.envelope->row_upper[199]=saved;mpq_clear(q);lab_model_free(&m);puts("source getters, masks and overlay identity controls passed");return 0;
}
