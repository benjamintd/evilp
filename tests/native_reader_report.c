#include "model.h"
#include "exact.h"
#include "common.h"
#include <stdio.h>
#include <sys/resource.h>
int main(int argc,char **argv){if(argc!=2)return 64;LabModel m={0};double start=lab_now();int bad=lab_model_read(&m,argv[1]);struct rusage r;getrusage(RUSAGE_SELF,&r);size_t zeros=0;if(m.native_exact)for(size_t k=0;k<m.nnz;k++)if(m.value[k]==0&&mpq_sgn(m.native_exact->value[k]))zeros++;printf("{\"read_ok\":%s,\"global_exact\":%s,\"native_exact\":%s,\"unsupported\":%s,\"rows\":%zu,\"cols\":%zu,\"nnz\":%zu,\"exact_nonzero_binary_zeros\":%zu,\"read_s\":%.9g,\"peak_rss_bytes\":%ld,\"error\":",bad?"false":"true",m.exact?"true":"false",m.native_exact?"true":"false",m.unsupported?"true":"false",m.rows,m.cols,m.nnz,zeros,lab_now()-start,r.ru_maxrss);lab_json_string(stdout,m.error);puts("}");lab_model_free(&m);return 0;}
