#include "model.h"
#include "exact.h"
#include "augment.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CHECK(x) do { if(!(x)){fprintf(stderr,"source scalar control %d: %s\n",__LINE__,#x);return 1;} } while(0)
int main(int argc,char **argv) {
 CHECK(argc==5);LabModel m={0},copy={0};mpq_t q,sentinel;mpq_inits(q,sentinel,NULL);mpq_set_ui(sentinel,7,9);
 CHECK(!lab_model_read(&m,argv[1]));CHECK(!m.exact&&m.native_exact&&m.objective_exact);
 double *x=calloc(m.cols,sizeof *x),*zero=calloc(m.cols,sizeof *zero),*a=calloc(m.cols,sizeof *a);CHECK(x&&zero&&a);x[0]=x[1]=1e15;
 CHECK(!lab_source_exact_objective(&m,0,q));CHECK(!lab_exact_decimal(sentinel,"1.000000000000000001")&&!mpq_cmp(q,sentinel));
 CHECK(!lab_source_exact_offset(&m,q)&&mpq_cmp_ui(q,1,10000)==0);
 CHECK(lab_source_objective_needs_exact(&m));double objective,violation,difference;
 CHECK(!lab_source_objective_difference(&m,x,zero,&difference)&&difference==.001);
 CHECK(!lab_source_objective_value(&m,x,&objective)&&objective==.0011);
 CHECK(!lab_exact_evaluate(&m,x,&objective,&violation)&&objective==.0011&&violation==0);
 LabAddedRow extra={"added",a,-INFINITY,1};LabAugmentLimits limits={4096,100000,1};
 CHECK(lab_model_append_rows(&m,1,&extra,&limits,&copy)==LAB_AUGMENT_OK);
 CHECK(copy.objective_exact&&copy.objective_exact!=m.objective_exact&&copy.envelope->objective!=m.envelope->objective);
 lab_model_free(&m);CHECK(!lab_exact_evaluate(&copy,x,&objective,&violation)&&objective==.0011&&violation==0);
 mpq_set_ui(q,7,9);copy.objective[0]=2;CHECK(lab_source_exact_objective(&copy,0,q)&&mpq_cmp_ui(q,7,9)==0);copy.objective[0]=1;
 lab_source_objective_free(copy.objective_exact);copy.objective_exact=NULL;mpq_set_ui(q,7,9);
 CHECK(lab_source_exact_objective(&copy,0,q)&&mpq_cmp_ui(q,7,9)==0);CHECK(lab_exact_evaluate(&copy,x,&objective,&violation)&&isnan(objective)&&isinf(violation));
 lab_model_free(&copy);free(x);free(zero);free(a);
 CHECK(!lab_model_read(&m,argv[2]));x=calloc(m.cols,sizeof *x);CHECK(x);x[0]=1e15;
 CHECK(m.native_exact&&!m.exact&&m.value[0]>0);CHECK(!lab_exact_evaluate(&m,x,&objective,&violation)&&violation>=.001&&violation<.001000000001);
 /* Restored native support must remain authoritative for an old binary-zero representation. */
 m.value[0]=0;CHECK(!lab_exact_evaluate(&m,x,&objective,&violation)&&violation>=.001&&violation<.001000000001);
 lab_exact_model_free(m.native_exact);m.native_exact=NULL;CHECK(!lab_exact_evaluate(&m,x,&objective,&violation)&&violation>=.001);
 m.row_lower[0]=-INFINITY;m.row_upper[0]=INFINITY;CHECK(!lab_exact_evaluate(&m,x,&objective,&violation)&&violation==0);
 lab_model_free(&m);free(x);
 CHECK(!lab_model_read(&m,argv[3]));x=calloc(m.cols,sizeof *x);CHECK(x);x[0]=1;
 CHECK(!lab_exact_evaluate(&m,x,&objective,&violation)&&violation>=1e-18&&violation<1e-14);
 m.col_lower[0]=.5;CHECK(!lab_exact_evaluate(&m,x,&objective,&violation)&&violation==0);lab_model_free(&m);free(x);
 CHECK(!lab_model_read(&m,argv[4]));CHECK(!m.exact&&!m.native_exact&&!m.objective_exact);x=calloc(m.cols,sizeof *x);CHECK(x);x[0]=1;
 mpq_set_ui(q,7,9);CHECK(lab_source_exact_objective(&m,0,q)&&mpq_cmp_ui(q,7,9)==0);
 CHECK(lab_exact_evaluate(&m,x,&objective,&violation));lab_model_free(&m);free(x);mpq_clears(q,sentinel,NULL);
 puts("source scalar, ambiguity, ownership and resource controls passed");return 0;
}
