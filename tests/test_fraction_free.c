#include "rational.h"
#include "common.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define CHECK(x) do {if(!(x)){fprintf(stderr,"exact solve contract line%d\n",__LINE__);exit(1);}}while(0)
static uint32_t state=0x21973a41;
static uint32_t random_value(void){state=1664525U*state+1013904223U;return state;}
static void generated(size_t n) {
 mpq_t *a=lab_alloc(n*n,sizeof(mpq_t)),*b=lab_alloc(n,sizeof(mpq_t)),*x=lab_alloc(n,sizeof(mpq_t));
 for(size_t k=0;k<n*n;k++)mpq_init(a[k]);
 for(size_t i=0;i<n;i++){mpq_init(b[i]);mpq_init(x[i]);mpq_set_si(x[i],(long)i-7,1+random_value()%17);mpq_canonicalize(x[i]);}
 mpq_t q;mpq_init(q);
 for(size_t i=0;i<n;i++){
  unsigned long denominator=1+random_value()%23;
  for(size_t j=0;j<n;j++){mpq_set_si(a[i*n+j],i==j?(long)(4*n+7):(long)(random_value()%7)-3,denominator);mpq_canonicalize(a[i*n+j]);}
 }
 /* Exercise pivot permutations without changing the exact known answer. */
 for(size_t i=0;i<n;i++){
  size_t row=random_value()%n;
  for(size_t j=0;j<n;j++)mpq_swap(a[i*n+j],a[row*n+j]);
 }
 for(size_t i=0;i<n;i++)for(size_t j=0;j<n;j++){mpq_mul(q,a[i*n+j],x[j]);mpq_add(b[i],b[i],q);}
 CHECK(lab_rational_solve(n,a,b,lab_now()+5)==0);
 for(size_t i=0;i<n;i++)CHECK(mpq_cmp(b[i],x[i])==0);
 mpq_clear(q);
 for(size_t k=0;k<n*n;k++)mpq_clear(a[k]);
 for(size_t i=0;i<n;i++){mpq_clear(b[i]);mpq_clear(x[i]);}free(a);free(b);free(x);
}
int main(void){
 for(size_t i=0;i<160;i++)generated(1+i%20);
 generated(64);generated(128);
 mpq_t a[4],b[2];for(size_t i=0;i<4;i++)mpq_init(a[i]);for(size_t i=0;i<2;i++)mpq_init(b[i]);
 mpq_set_ui(a[1],1,1);mpq_set_ui(a[2],2,1);mpq_set_ui(a[3],3,1);mpq_set_ui(b[0],2,1);mpq_set_ui(b[1],8,1);
 CHECK(lab_rational_solve(2,a,b,lab_now()+1)==0);CHECK(mpq_cmp_ui(b[0],1,1)==0&&mpq_cmp_ui(b[1],2,1)==0);
 for(size_t i=0;i<4;i++)mpq_set_ui(a[i],0,1);
 CHECK(lab_rational_solve(2,a,b,lab_now()+1)!=0);
 CHECK(lab_rational_solve(2,a,b,lab_now()-1)!=0);
 CHECK(lab_rational_solve(2,a,b,NAN)!=0);
 CHECK(lab_rational_solve(129,NULL,NULL,lab_now()+1)!=0);
 mpz_setbit(mpq_numref(a[0]),LAB_RATIONAL_MAX_BITS+1);
 CHECK(lab_rational_solve(2,a,b,lab_now()+1)!=0);
 for(size_t i=0;i<4;i++)mpq_clear(a[i]);for(size_t i=0;i<2;i++)mpq_clear(b[i]);
 puts("162 exact known-solution systems and pivot/cap/singular controls passed");return 0;
}
