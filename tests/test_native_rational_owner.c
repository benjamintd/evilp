#include "model.h"
#include "exact.h"
#include "augment.h"
#include "root_cut.h"
#include "common.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"native owner control%d: %s\n",__LINE__,#c);return 1;}}while(0)
int main(int argc,char **argv){
 CHECK(argc==3);LabModel source={0},owned={0};CHECK(!lab_model_read(&source,argv[1]));CHECK(!source.exact&&source.native_exact&&source.rows==200&&source.cols==300);
 size_t zero=0;for(size_t k=0;k<source.nnz;k++)if(source.value[k]==0&&mpq_sgn(source.native_exact->value[k]))zero++;
 double *a=calloc(source.cols,sizeof *a),*lo=malloc(source.cols*sizeof *lo),*hi=malloc(source.cols*sizeof *hi),*rl=malloc(source.rows*sizeof *rl),*ru=malloc(source.rows*sizeof *ru);CHECK(a&&lo&&hi&&rl&&ru);
 memcpy(lo,source.col_lower,source.cols*sizeof *lo);memcpy(hi,source.col_upper,source.cols*sizeof *hi);lo[0]=.25;
 memcpy(rl,source.row_lower,source.rows*sizeof *rl);memcpy(ru,source.row_upper,source.rows*sizeof *ru);rl[199]=-INFINITY;ru[199]=INFINITY;
 LabModel view=source;view.col_lower=lo;view.col_upper=hi;view.row_lower=rl;view.row_upper=ru;
 LabAddedRow added={"extra_constant",a,-1,INFINITY};LabAugmentLimits cap={4096,100000,4};CHECK(lab_model_append_rows(&view,1,&added,&cap,&owned)==LAB_AUGMENT_OK);CHECK(owned.native_exact&&owned.native_exact!=source.native_exact&&!owned.exact);
 CHECK(owned.native_exact->col_lower_source[0]==source.native_exact->col_lower_source[0]);for(size_t k=0;k<source.nnz;k++)CHECK(!mpq_cmp(source.native_exact->value[k],owned.native_exact->value[k]));
 mpq_t q,original;mpq_inits(q,original,NULL);mpq_set(original,source.native_exact->col_lower[0]);CHECK(!lab_source_exact_column_bound(&owned,0,0,q)&&mpq_cmp_ui(q,1,4)==0);CHECK(!lab_source_exact_column_value(&owned,0,0,source.col_lower[0],q)&&!mpq_cmp(q,original));
 lab_model_free(&source);free(a);free(lo);free(hi);free(rl);free(ru);
 LabLpOptions o={.time_limit=30,.iteration_limit=10000,.max_dimension=4096};LabLpSession *session=lab_lp_session_create(&owned,&o);CHECK(session);LabLpResult r;lab_lp_session_solve(session,&r);CHECK(!strcmp(r.status,"optimal")&&r.has_primal);CHECK(r.exact_primal_attempts==0);
 LabRootCutBatch b;lab_root_cut_batch_init(&b);LabRootCutLimits limits=lab_root_cut_limits(lab_now()+10);CHECK(lab_root_cut_derive(session,&owned,&limits,&b)==LAB_ROOT_CUT_OK&&b.count==1);
 FILE *f=fopen(argv[2],"wx");CHECK(f&&fputs(b.cut[0].certificate,f)>=0&&!fclose(f));lab_root_cut_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(session);lab_model_free(&owned);mpq_clears(q,original,NULL);printf("{\"passed\":true,\"native_zero_coefficients\":%zu,\"global_exact_primal_attempts\":0}\n",zero);return 0;
}
