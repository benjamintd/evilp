#include "root_cut.h"
#include "common.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(c) do { if (!(c)) { fprintf(stderr,"root-cut contract line %d: %s\n",__LINE__,#c); return 1; } } while (0)

/* The Python oracle creates fresh source/proof paths and independently replays
 * this producer's certificates with milpcheck. Optimization remains in C. */
int main(int argc, char **argv) {
    if (argc < 3 || argc > 8) return 64;
    LabModel model = {0};
    CHECK(!lab_model_read(&model, argv[1]));
    LabLpOptions options = {.time_limit=30,.iteration_limit=10000,.max_dimension=4096};
    if (argc >= 4) options.dense = atoi(argv[3]);
    if (argc >= 5) options.no_scaling = atoi(argv[4]);
    if (argc >= 8) options.no_exact = atoi(argv[7]);
    LabLpSession *session = lab_lp_session_create(&model, &options);
    CHECK(session);
    LabRootCutLimits limits = lab_root_cut_limits(lab_now()+20);
    LabRootCutBatch batch;
    lab_root_cut_batch_init(&batch);
    CHECK(lab_root_cut_derive(session,&model,&limits,&batch)==LAB_ROOT_CUT_STALE);
    if (argc >= 7) CHECK(!lab_lp_session_set_bounds(session,0,strtod(argv[5],NULL),strtod(argv[6],NULL)));
    LabLpResult result;
    lab_lp_session_solve(session,&result);
    CHECK(!strcmp(result.status,"optimal")&&result.has_primal);
    LabRootCutStatus status=lab_root_cut_derive(session,&model,&limits,&batch);
    CHECK(status==LAB_ROOT_CUT_OK||status==LAB_ROOT_CUT_NONE||status==LAB_ROOT_CUT_SOURCE_UNAVAILABLE);
    size_t controls=1;
    if (status==LAB_ROOT_CUT_OK) {
        CHECK(batch.count&&batch.count<=4&&batch.columns==model.cols&&batch.box);
        LabRootCutBatch saved=batch;
        LabRootCutLimits bad=limits;
        bad.native.deadline=lab_now()-1;
        CHECK(lab_root_cut_derive(session,&model,&bad,&batch)==LAB_ROOT_CUT_DEADLINE);
        CHECK(!memcmp(&saved,&batch,sizeof batch));controls++;
        bad=limits;bad.max_certificate_bytes=1;
        CHECK(lab_root_cut_derive(session,&model,&bad,&batch)==LAB_ROOT_CUT_RESOURCE);
        CHECK(!memcmp(&saved,&batch,sizeof batch));controls++;
        bad=limits;bad.max_certificate_bytes=0;
        CHECK(lab_root_cut_derive(session,&model,&bad,&batch)==LAB_ROOT_CUT_INVALID);
        CHECK(!memcmp(&saved,&batch,sizeof batch));controls++;
        unsigned char *integer=calloc(model.cols?model.cols:1,1);CHECK(integer);
        memcpy(integer,model.integer,model.cols);integer[0]=(unsigned char)!integer[0];
        LabModel other=model;other.integer=integer;
        CHECK(lab_root_cut_derive(session,&other,&limits,&batch)==LAB_ROOT_CUT_SOURCE_UNAVAILABLE);
        CHECK(!memcmp(&saved,&batch,sizeof batch));free(integer);controls++;
        if (model.nnz) {
            double *values=malloc(model.nnz*sizeof *values);CHECK(values);
            memcpy(values,model.value,model.nnz*sizeof *values);
            values[0]=nextafter(values[0],INFINITY);other=model;
            other.value=values;other.exact=NULL;other.parsed_decimal=0;
            CHECK(lab_root_cut_derive(session,&other,&limits,&batch)==LAB_ROOT_CUT_SOURCE_UNAVAILABLE);
            CHECK(!memcmp(&saved,&batch,sizeof batch));free(values);controls++;
        }
        double lo=argc>=7?strtod(argv[5],NULL):model.col_lower[0];
        double hi=argc>=7?strtod(argv[6],NULL):model.col_upper[0];
        CHECK(!lab_lp_session_set_bounds(session,0,lo,hi));
        CHECK(lab_root_cut_derive(session,&model,&limits,&batch)==LAB_ROOT_CUT_OK);controls++;
        if (isfinite(lo)&&isfinite(hi)&&lo<hi) {
            CHECK(!lab_lp_session_set_bounds(session,0,nextafter(lo,hi),hi));saved=batch;
            CHECK(lab_root_cut_derive(session,&model,&limits,&batch)==LAB_ROOT_CUT_STALE);
            CHECK(!memcmp(&saved,&batch,sizeof batch));controls++;
        }
    }
    size_t columns=model.cols;
    /* The batch owns its proof, coefficient and box allocations. */
    lab_lp_result_free(&result);lab_lp_session_free(session);lab_model_free(&model);
    printf("{\"status\":");lab_json_string(stdout,lab_root_cut_status_name(status));
    printf(",\"cuts\":%zu,\"controls\":%zu,\"columns\":%zu,\"published\":[",batch.count,controls,columns);
    mpq_t q;mpq_init(q);
    for (size_t i=0;i<batch.count;++i) {
        char path[4096];int length=snprintf(path,sizeof path,"%s-%zu.cert",argv[2],i);
        CHECK(length>0&&(size_t)length<sizeof path);FILE *f=fopen(path,"wx");CHECK(f);
        CHECK(fputs(batch.cut[i].certificate,f)>=0);CHECK(!fclose(f));
        CHECK(batch.cut[i].columns==columns&&isfinite(batch.cut[i].lower));
        printf("%s{\"certificate\":",i?",":"");lab_json_string(stdout,path);fputs(",\"lower\":",stdout);
        mpq_set_d(q,batch.cut[i].lower);char *text=mpq_get_str(NULL,10,q);lab_json_string(stdout,text);free(text);
        fputs(",\"coefficient\":[",stdout);
        for (size_t j=0;j<columns;++j) {
            CHECK(isfinite(batch.cut[i].coefficient[j]));if(j)fputc(',',stdout);
            mpq_set_d(q,batch.cut[i].coefficient[j]);text=mpq_get_str(NULL,10,q);lab_json_string(stdout,text);free(text);
        }
        fputs("]}",stdout);
    }
    mpq_clear(q);puts("]}");lab_root_cut_batch_free(&batch);lab_root_cut_batch_free(&batch);
    return 0;
}
