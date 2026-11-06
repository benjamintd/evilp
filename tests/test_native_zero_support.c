#include "model.h"
#include "exact.h"
#include "native_tableau.h"
#include "root_cut.h"
#include "common.h"
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#define CHECK(c) do{if(!(c)){fprintf(stderr,"zero support%d: %s\n",__LINE__,#c);return 1;}}while(0)
int main(int argc,char **argv){
 CHECK(argc==2);LabModel m={0};CHECK(!lab_model_read(&m,argv[1]));CHECK(!m.exact&&m.native_exact);
 size_t source=SIZE_MAX;for(size_t k=m.col_start[0];k<m.col_start[1];k++)if(m.row_index[k]==199&&mpq_sgn(m.native_exact->value[k]))source=k;CHECK(source!=SIZE_MAX&&m.value[source]!=0);
 /* The reader repairs this entry. Simulate a preexisting binary-zero model to
  * retain the exporter regression for restoring exact canonical/basis support. */
 m.value[source]=0;
 for(int scaling=0;scaling<2;scaling++) {
  LabLpOptions o={.time_limit=30,.iteration_limit=10000,.max_dimension=4096,.no_scaling=!scaling};LabLpSession*s=lab_lp_session_create(&m,&o);CHECK(s);LabLpResult r;lab_lp_session_solve(s,&r);CHECK(!strcmp(r.status,"optimal"));
  LabNativeBatch b;lab_native_batch_init(&b);LabRootCutLimits limits=lab_root_cut_limits(lab_now()+10);CHECK(lab_lp_session_export_native(s,&limits.native,&b)==LAB_NATIVE_OK);
  size_t row=SIZE_MAX,pivot=SIZE_MAX;for(size_t i=0;i<b.dimension;i++){if(b.row_id[i]==m.row_index[source])row=i;if(b.basis[i]==0)pivot=i;}CHECK(row!=SIZE_MAX&&pivot!=SIZE_MAX);
  int canonical=0,basis=0;mpq_t expected;mpq_init(expected);mpq_mul(expected,m.native_exact->value[source],b.tableau[0].signed_basic_scale);
  for(size_t k=b.canonical_matrix.start[0];k<b.canonical_matrix.start[1];k++)if(b.canonical_matrix.row[k]==row){CHECK(!mpq_cmp(expected,b.canonical_matrix.value[k]));canonical++;}
  for(size_t k=b.basis_matrix.start[pivot];k<b.basis_matrix.start[pivot+1];k++)if(b.basis_matrix.row[k]==row){CHECK(!mpq_cmp(expected,b.basis_matrix.value[k]));basis++;}
  CHECK(canonical==1&&basis==1);mpq_clear(expected);lab_native_batch_free(&b);lab_lp_result_free(&r);lab_lp_session_free(s);
 }
 lab_model_free(&m);puts("exact support retained in canonical and basis matrices for both scales");return 0;
}
