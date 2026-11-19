/* The private bound is checked against exact GMP products and the retained
 * compensated calculation, including cancellation and extreme magnitudes. */
#include "../src/linear.c"
#include <stdio.h>
#include <gmp.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"residual bound line%d\n",__LINE__);return 1;}}while(0)
static uint64_t rng=UINT64_C(0xa75fc3d12920e8b7);
static uint64_t word(void){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;return rng;}
static int verify(LabLinear *s, const double *rhs, const double *x, int transpose, int *accepted){
 size_t n=s->n;double *r=lab_alloc(n,sizeof(double));double error;
 if(!fast_residual_bound(s,rhs,x,r,transpose,&error)){free(r);return 0;}
 (*accepted)++;CHECK(error<=1e-14);
 double full;CHECK(residual(s,rhs,x,r,transpose,&full)==0);CHECK(full<=error);
 mpq_t *sum=lab_alloc(n,sizeof(mpq_t)),*den=lab_alloc(n,sizeof(mpq_t));
 mpq_t a,b,t,q,bound;mpq_inits(a,b,t,q,bound,NULL);mpq_set_d(bound,error);
 for(size_t i=0;i<n;i++){mpq_init(sum[i]);mpq_init(den[i]);mpq_set_d(sum[i],rhs[i]);mpq_abs(den[i],sum[i]);mpq_set_ui(t,1,1);mpq_add(den[i],den[i],t);}
 for(size_t j=0;j<n;j++){
  size_t count;const size_t *row;const double *v;current_column(s,j,&count,&row,&v);
  for(size_t k=0;k<count;k++){size_t i=row[k],at=transpose?j:i;mpq_set_d(a,v[k]);mpq_set_d(b,x[transpose?i:j]);mpq_mul(t,a,b);mpq_sub(sum[at],sum[at],t);mpq_abs(t,t);mpq_add(den[at],den[at],t);}
 }
 for(size_t i=0;i<n;i++){mpq_abs(q,sum[i]);mpq_div(q,q,den[i]);CHECK(mpq_cmp(q,bound)<=0);mpq_clear(sum[i]);mpq_clear(den[i]);}
 free(sum);free(den);mpq_clears(a,b,t,q,bound,NULL);free(r);return 0;
}
int main(void){
 int accepted=0;
 for(size_t trial=0;trial<2000;trial++){
  size_t n=1+word()%24,nnz=0;size_t *start=lab_alloc(n+1,sizeof(size_t)),*row=lab_alloc(n*n,sizeof(size_t));double *value=lab_alloc(n*n,sizeof(double)),*x=lab_alloc(n,sizeof(double)),*b=lab_alloc(n,sizeof(double));
  for(size_t i=0;i<n;i++)x[i]=scalbn((double)((int)(word()%31)-15),(int)(word()%40)-20);
  for(size_t j=0;j<n;j++){start[j]=nnz;for(size_t i=0;i<n;i++)if(i==j||word()%3==0){row[nnz]=i;value[nnz++]=scalbn((double)((int)(word()%31)-15),(int)(word()%40)-20);}}start[n]=nnz;
  LabLinear s={.n=n,.base_col_start=start,.base_row_index=row,.base_value=value};
  for(int transpose=0;transpose<2;transpose++){
   memset(b,0,n*sizeof(double));for(size_t j=0;j<n;j++)for(size_t k=start[j];k<start[j+1];k++){size_t i=row[k];b[transpose?j:i]=fma(value[k],x[transpose?i:j],b[transpose?j:i]);}
   if(trial%7==0)b[0]+=1e-4;
   CHECK(verify(&s,b,x,transpose,&accepted)==0);
  }
  free(start);free(row);free(value);free(x);free(b);
 }
 CHECK(accepted>2500);
 size_t start[]={0,1},row[]={0};double a[]={DBL_MIN},x[]={DBL_MIN},rhs[]={0},r[1],e;LabLinear s={.n=1,.base_col_start=start,.base_row_index=row,.base_value=a};
 CHECK(!fast_residual_bound(&s,rhs,x,r,0,&e));a[0]=DBL_MAX;x[0]=2;CHECK(!fast_residual_bound(&s,rhs,x,r,0,&e));a[0]=1;x[0]=NAN;CHECK(!fast_residual_bound(&s,rhs,x,r,0,&e));
 a[0]=.1;x[0]=.3;rhs[0]=.03;CHECK(verify(&s,rhs,x,0,&accepted)==0);
 CHECK(fesetround(FE_UPWARD)==0);CHECK(!fast_residual_bound(&s,rhs,x,r,0,&e));CHECK(fesetround(FE_TONEAREST)==0);
 printf("4000 exact residual comparisons; %d accepted bounds; underflow/overflow/rounding controls passed\n",accepted);return 0;
}
