#include "pivot.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
static void require(int condition, const char *why) {
    if (!condition) { fprintf(stderr, "scaled ratio: %s\n", why); exit(1); }
}
int main(void) {
    size_t basis[2] = {1,2}, row;
    size_t cases = 0, old_unsafe = 0;
    /* Exact dyadic limiting steps differ by 2^-56. A dimensionless 1e-12
     * tie incorrectly chooses row0 and crosses row1's original-unit bound.
     * Covary both row units independently; the true permissible step stays
     * unchanged and each original error allowance transforms with the row. */
    for (int first = -8; first <= 8; first += 4)
        for (int second = -8; second <= 8; second += 4) {
            double x[2] = {scalbn(0x1p-4,first),scalbn(0x1p-38-0x1p-60,second)};
            double d[2] = {scalbn(0x1p30,first),scalbn(0x1p-4,second)};
            double t[2] = {scalbn(0x1p-40,first),scalbn(0x1p-70,second)};
            for (int bland = 0; bland <= 1; ++bland) {
                require(!lab_primal_ratio_scaled(2,x,d,basis,t,bland,&row) && row==1,
                        "unit changes preserve the genuinely limiting row");
                ++cases;
            }
            require(!lab_primal_ratio(2,x,d,basis,1,&row), "legacy comparison row");
            if (row == 0 && fmal(-((long double)x[row]/d[row]),d[1],x[1]) < -t[1])
                ++old_unsafe;
        }
    require(old_unsafe==25,"legacy absolute ties hide material original violations");
    double x[2]={0,1},d[2]={1e-13,1},t[2]={0,1e-8};
    require(lab_primal_ratio_scaled(2,x,d,basis,t,1,&row)==-1,
            "a too-weak positive direction remains a limiting constraint");
    d[0]=d[1]=0;
    require(lab_primal_ratio_scaled(2,x,d,basis,t,0,&row)==1,"truly unbounded direction");
    d[0]=d[1]=1;x[0]=x[1]=1;t[0]=t[1]=1e-10;basis[0]=9;basis[1]=2;
    for(int bland=0;bland<=1;++bland)
        require(!lab_primal_ratio_scaled(2,x,d,basis,t,bland,&row)&&row==1,
                "exact ratio ties remain deterministic");
    t[0]=NAN;
    require(lab_primal_ratio_scaled(2,x,d,basis,t,0,&row)==-1,"nonfinite allowance rejected");
    t[0]=-1;
    require(lab_primal_ratio_scaled(2,x,d,basis,t,0,&row)==-1,"negative allowance rejected");
    t[0]=2*LAB_PRIMAL_TOLERANCE;
    require(lab_primal_ratio_scaled(2,x,d,basis,t,0,&row)==-1,"no relaxed global tolerance");
    printf("{\"scaled_ratio_cases\":%zu,\"legacy_unsafe_cases\":%zu}\n",cases,old_unsafe);
    return 0;
}
