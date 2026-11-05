#define _DEFAULT_SOURCE
#define _DARWIN_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#include "common.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <signal.h>
#include <fcntl.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>

typedef struct {char *id,*path,*sha,*reference,*value;} Instance;
typedef struct {char *name,*binary,*sha,*version;} Solver;
typedef struct {int code,signal,timed_out;double wall,cpu;long rss;} Process;
static void die(const char *s){fprintf(stderr,"milpbench: %s\n",s);exit(2);}
static char *pathcat(const char *a,const char *b){size_t n=strlen(a)+strlen(b)+2;char *p=lab_alloc(n,1);snprintf(p,n,"%s/%s",a,b);return p;}
static char *absolute(const char *s){char *p=realpath(s,NULL);if(!p)die(s);return p;}
static void mkdirs(const char *path){char *s=lab_strdup(path);for(char *p=s+1;*p;p++)if(*p=='/'){*p=0;if(mkdir(s,0755)&&errno!=EEXIST)die("mkdir");*p='/';}if(mkdir(s,0755)&&errno!=EEXIST)die("mkdir");free(s);}
static int safe_id(const char *s){if(!*s||!strcmp(s,".")||!strcmp(s,".."))return 0;for(;*s;s++)if(!isalnum((unsigned char)*s)&&!strchr("._-",*s))return 0;return 1;}
static void verify(const char *path,const char *expected){char h[65];if(strlen(expected)!=64||lab_sha256(path,h)||strcmp(h,expected))die("SHA-256 mismatch or unreadable locked file");}
static Process run(char *const args[],const char *log,double limit){
    Process r={0};double begin=lab_now();pid_t pid=fork();if(pid<0)die("fork");
    if(!pid){
        if(setpgid(0,0))_exit(126);
        int fd=open(log,O_WRONLY|O_CREAT|O_EXCL,0644);if(fd<0)_exit(126);
        if(dup2(fd,STDOUT_FILENO)<0||dup2(fd,STDERR_FILENO)<0)_exit(126);close(fd);
        int input=open("/dev/null",O_RDONLY);if(input>=0){dup2(input,STDIN_FILENO);close(input);}
        setenv("OMP_NUM_THREADS","1",1);setenv("OPENBLAS_NUM_THREADS","1",1);setenv("MKL_NUM_THREADS","1",1);
        execv(args[0],args);perror("execv");_exit(127);
    }
    (void)setpgid(pid,pid);int status=0;struct rusage usage;memset(&usage,0,sizeof usage);
    int terminated=0;double kill_at=0;
    for(;;){
        pid_t p=wait4(pid,&status,WNOHANG,&usage);if(p==pid)break;if(p<0&&errno!=EINTR)die("wait4");
        double elapsed=lab_now()-begin;
        if(!terminated&&elapsed>limit){r.timed_out=1;terminated=1;kill_at=elapsed+1;kill(-pid,SIGTERM);}
        if(terminated==1&&elapsed>kill_at){kill(-pid,SIGKILL);terminated=2;}
        struct timespec pause={0,20000000};nanosleep(&pause,NULL);
    }
    /* Also clean up any descendants that outlive the main process. */
    kill(-pid,SIGKILL);
    r.wall=lab_now()-begin;r.code=WIFEXITED(status)?WEXITSTATUS(status):-1;r.signal=WIFSIGNALED(status)?WTERMSIG(status):0;
    r.cpu=(double)usage.ru_utime.tv_sec+(double)usage.ru_utime.tv_usec/1e6+(double)usage.ru_stime.tv_sec+(double)usage.ru_stime.tv_usec/1e6;
#ifdef __APPLE__
    r.rss=usage.ru_maxrss/1024;
#else
    r.rss=usage.ru_maxrss;
#endif
    return r;
}
static char *readfile(const char *path){FILE *f=fopen(path,"rb");if(!f)return lab_strdup("");if(fseek(f,0,SEEK_END)){fclose(f);return lab_strdup("");}long n=ftell(f);if(n<0||n>16*1024*1024)die("oversized result");rewind(f);char *s=lab_alloc((size_t)n+1,1);size_t got=fread(s,1,(size_t)n,f);s[got]=0;fclose(f);return s;}
static void argv_json(FILE *f,char *const args[]){fputc('[',f);for(int i=0;args[i];i++){if(i)fputc(',',f);lab_json_string(f,args[i]);}fputc(']',f);}
static const char *status_from_log(const char *solver,const char *path,Process p,char bound[128]){
    FILE *f=fopen(path,"r");if(!f)return "launch_error";char line[16384];const char *status="unknown";bound[0]=0;int highs_report=0;
    while(fgets(line,sizeof line,f)){
        if(!strcmp(solver,"highs")&&strstr(line,"Solving report"))highs_report=1;
        if(!strcmp(solver,"highs")&&(strstr(line,"Model status")||(highs_report&&!strncmp(line,"  Status ",9)))){
            if(strstr(line,"Optimal"))status="optimal";else if(strstr(line,"Infeasible"))status="infeasible";else if(strstr(line,"Unbounded"))status="unbounded";else if(strstr(line,"Time limit"))status="time_limit";
        }else if(!strcmp(solver,"scip")&&strstr(line,"SCIP Status")){
            if(strstr(line,"optimal solution"))status="optimal";else if(strstr(line,"infeasible"))status="infeasible";else if(strstr(line,"unbounded"))status="unbounded";else if(strstr(line,"time limit"))status="time_limit";
        }else if(!strcmp(solver,"cbc")){
            if(strstr(line,"Result - Optimal solution found"))status="optimal";
            else if(strstr(line,"Result - Problem proven infeasible")||strstr(line,"Problem is infeasible"))status="infeasible";
            else if(strstr(line,"Result - Linear relaxation unbounded"))status="unbounded";
            else if(strstr(line,"Result - Stopped on time limit"))status="time_limit";
        }
        char *start=NULL;
        if(!strcmp(solver,"highs"))start=strstr(line,"Dual bound");
        else if(!strcmp(solver,"scip"))start=strstr(line,"Dual Bound");
        else start=strstr(line,"Lower bound:");
        if(start){while(*start&&!isdigit((unsigned char)*start)&&*start!='-'&&*start!='+')start++;char *end;double d=strtod(start,&end);if(end!=start&&isfinite(d))snprintf(bound,128,"%.17g",d);}
    }
    fclose(f);
    if(p.timed_out)return "watchdog_timeout";if(p.code==127||p.code==126)return "launch_error";if(p.signal||p.code!=0)return "process_error";return status;
}
static int split(char *s,char **v,int max){size_t n=strlen(s);while(n&&(s[n-1]=='\n'||s[n-1]=='\r'))s[--n]=0;int count=0;char *p=s;for(;;){if(count==max)return -1;v[count++]=p;char *t=strchr(p,'\t');if(!t)break;*t=0;p=t+1;}return count;}
static size_t load_instances(const char *path,Instance **out){
    FILE *f=fopen(path,"r");if(!f)die("cannot open suite");char *line=NULL,*v[5];size_t size=0,n=0;ssize_t len;
    Instance *a=NULL;
    while((len=getline(&line,&size,f))>=0){if(len<=1||line[0]=='#'||!strncmp(line,"id\t",3))continue;if(split(line,v,5)!=5||!safe_id(v[0]))die("malformed suite");
        Instance *next=realloc(a,(n+1)*sizeof *a);if(!next)die("allocation");a=next;
        for(size_t j=0;j<n;j++)if(!strcmp(a[j].id,v[0]))die("duplicate instance ID");
        a[n]=(Instance){lab_strdup(v[0]),absolute(v[1]),lab_strdup(v[2]),lab_strdup(v[3]),lab_strdup(v[4])};n++;
    }free(line);fclose(f);*out=a;return n;
}
static size_t load_solvers(const char *path,Solver **out){
    FILE *f=fopen(path,"r");if(!f)die("cannot open solvers config");char *line=NULL,*v[4];size_t size=0,n=0;Solver *a=NULL;
    while(getline(&line,&size,f)>=0){if(line[0]=='#'||!strncmp(line,"solver\t",7)||line[0]=='\n')continue;if(split(line,v,4)!=4)die("malformed solvers config");
        if(strcmp(v[0],"highs")&&strcmp(v[0],"scip")&&strcmp(v[0],"cbc"))die("unsupported solver adapter");
        Solver *next=realloc(a,(n+1)*sizeof *a);if(!next)die("allocation");a=next;
        for(size_t j=0;j<n;j++)if(!strcmp(a[j].name,v[0]))die("duplicate solver");
        a[n]=(Solver){lab_strdup(v[0]),absolute(v[1]),lab_strdup(v[2]),lab_strdup(v[3])};verify(a[n].binary,a[n].sha);n++;
    }free(line);fclose(f);*out=a;return n;
}
static void one(const Instance *inst,const Solver *s,const char *checker,const char *checker_sha,const char *runner_sha,const char *out,double limit,double grace,int seed,FILE *all){
    char id[1024];snprintf(id,sizeof id,"%s/%s/seed-%d",inst->id,s->name,seed);char *dir=pathcat(out,id);mkdirs(dir);
    char *log=pathcat(dir,"solver.log"),*sol=pathcat(dir,"solution.sol"),*options=pathcat(dir,"settings.txt"),*result=pathcat(dir,"run.json"),*checklog=pathcat(dir,"checker.json");
    char secs[64],seedstr[32];snprintf(secs,sizeof secs,"%.9g",limit);snprintf(seedstr,sizeof seedstr,"%d",seed);
    char *args[40]={0};int k=0;args[k++]=s->binary;char commands[11][PATH_MAX+128];memset(commands,0,sizeof commands);
    FILE *opt=fopen(options,"wx");if(!opt)die("result directory is not new");
    if(!strcmp(s->name,"highs")){
        fprintf(opt,"time_limit = %s\nthreads = 1\nrandom_seed = %s\nmip_rel_gap = 0\nmip_abs_gap = 0\nwrite_solution_to_file = true\nwrite_solution_style = 0\nsolution_file = %s\nlog_file = %s/highs.log\n",secs,seedstr,sol,dir);
        args[k++]=inst->path;args[k++]="--options_file";args[k++]=options;
    }else if(!strcmp(s->name,"cbc")){
        args[k++]=inst->path;args[k++]="-threads";args[k++]="1";args[k++]="-randomSeed";args[k++]=seedstr;args[k++]="-randomCbcSeed";args[k++]=seedstr;
        args[k++]="-timeMode";args[k++]="elapsed";
        args[k++]="-seconds";args[k++]=secs;args[k++]="-ratioGap";args[k++]="0";args[k++]="-allowableGap";args[k++]="0";
        args[k++]="-solve";args[k++]="-solution";args[k++]=sol;args[k++]="-quit";
        fprintf(opt,"threads=1\nrandomSeed=%s\nrandomCbcSeed=%s\ntimeMode=elapsed\nseconds=%s\nratioGap=0\nallowableGap=0\n",seedstr,seedstr,secs);
    }else{
        if(strpbrk(inst->path,"\"\n\r")||strpbrk(sol,"\"\n\r"))die("SCIP adapter cannot represent quotes/newlines in filenames");
        snprintf(commands[0],sizeof commands[0],"set limits time %s",secs);
        snprintf(commands[1],sizeof commands[1],"set parallel maxnthreads 1");
        snprintf(commands[2],sizeof commands[2],"set randomization randomseedshift %s",seedstr);
        snprintf(commands[3],sizeof commands[3],"set limits gap 0");
        snprintf(commands[4],sizeof commands[4],"set limits absgap 0");
        snprintf(commands[5],sizeof commands[5],"read \"%s\"",inst->path);
        snprintf(commands[6],sizeof commands[6],"optimize");
        snprintf(commands[7],sizeof commands[7],"write solution \"%s\"",sol);
        snprintf(commands[8],sizeof commands[8],"quit");
        snprintf(commands[9],sizeof commands[9],"set timing clocktype 2");
        snprintf(commands[10],sizeof commands[10],"set presolving milp threads 1");
        args[k++]="-c";args[k++]=commands[9];fprintf(opt,"%s\n",commands[9]);
        args[k++]="-c";args[k++]=commands[10];fprintf(opt,"%s\n",commands[10]);
        for(int i=0;i<9;i++){args[k++]="-c";args[k++]=commands[i];fprintf(opt,"%s\n",commands[i]);}
    }
    fclose(opt);args[k]=NULL;
    Process p=run(args,log,limit+grace);char bound[128];const char *claimed=status_from_log(s->name,log,p,bound);
    Process cp={0};int has_solution=access(sol,R_OK)==0;char *check=NULL;
    if(has_solution){char *ca[]={(char*)checker,"--model",inst->path,"--solution",sol,"--format",s->name,NULL};cp=run(ca,checklog,120);check=readfile(checklog);}
    if(!check||check[0]!='{'||cp.timed_out||cp.signal||cp.code<0||cp.code>4){free(check);check=lab_strdup("{\"status\":\"unavailable\",\"checked_feasible\":false}");}
    size_t checklen=strlen(check);while(checklen&&isspace((unsigned char)check[checklen-1]))check[--checklen]=0;
    FILE *f=fopen(result,"wx");if(!f)die("cannot write run result");
    fprintf(f,"{\"schema\":1,\"instance\":");lab_json_string(f,inst->id);fprintf(f,",\"model_path\":");lab_json_string(f,inst->path);
    fprintf(f,",\"model_sha256\":");lab_json_string(f,inst->sha);fprintf(f,",\"solver\":");lab_json_string(f,s->name);
    fprintf(f,",\"solver_version\":");lab_json_string(f,s->version);fprintf(f,",\"solver_sha256\":");lab_json_string(f,s->sha);
    fprintf(f,",\"checker_sha256\":");lab_json_string(f,checker_sha);fprintf(f,",\"runner_sha256\":");lab_json_string(f,runner_sha);
    fprintf(f,",\"seed\":%d,\"requested_threads\":1,\"time_limit_s\":%.9g,\"watchdog_grace_s\":%.9g,\"wall_s\":%.9g,\"cpu_s\":%.9g,\"peak_rss_kib\":%ld,\"exit_code\":%d,\"signal\":%d,\"watchdog_timeout\":%s",seed,limit,grace,p.wall,p.cpu,p.rss,p.code,p.signal,p.timed_out?"true":"false");
    fprintf(f,",\"claimed_status\":");lab_json_string(f,claimed);fprintf(f,",\"reported_dual_bound\":");if(bound[0])lab_json_string(f,bound);else fputs("null",f);
    fprintf(f,",\"reference_status\":");lab_json_string(f,inst->reference);fprintf(f,",\"reference_objective\":");lab_json_string(f,inst->value);
    fprintf(f,",\"argv\":");argv_json(f,args);fprintf(f,",\"checker_exit_code\":");if(has_solution)fprintf(f,"%d",cp.code);else fputs("null",f);
    fprintf(f,",\"checker_wall_s\":%.9g,\"checker\":%s}\n",cp.wall,check);fclose(f);
    char *json=readfile(result);fputs(json,all);fflush(all);printf("%s %s seed=%d status=%s wall=%.3fs check=%s\n",inst->id,s->name,seed,claimed,p.wall,has_solution?(cp.code==0?"valid":cp.code==1?"invalid":"unavailable"):"no_solution");fflush(stdout);
    free(json);free(check);free(dir);free(log);free(sol);free(options);free(result);free(checklog);
}
int main(int argc,char **argv){
    const char *suite=NULL,*config="config/solvers.tsv",*output=NULL,*checker_arg="build/milpcheck";double limit=2,grace=15;int seed=1;size_t max=SIZE_MAX;
    for(int i=1;i<argc;i++){
        if(i+1>=argc)die("missing argument value");
        const char *key=argv[i],*value=argv[++i];char *end;
        if(!strcmp(key,"--suite"))suite=value;else if(!strcmp(key,"--solvers"))config=value;else if(!strcmp(key,"--out"))output=value;else if(!strcmp(key,"--checker"))checker_arg=value;
        else if(!strcmp(key,"--limit")){limit=strtod(value,&end);if(*end||!isfinite(limit)||limit<=0||limit>86400)die("invalid time limit");}
        else if(!strcmp(key,"--grace")){grace=strtod(value,&end);if(*end||!isfinite(grace)||grace<0||grace>3600)die("invalid grace");}
        else if(!strcmp(key,"--seed")){long x=strtol(value,&end,10);if(*end||x<1||x>INT_MAX)die("seed must be 1..INT_MAX; zero is wall-clock seeded in some solvers");seed=(int)x;}
        else if(!strcmp(key,"--max-instances")){unsigned long x=strtoul(value,&end,10);if(*end||!x)die("invalid instance count");max=(size_t)x;}
        else die("unknown option");
    }
    if(!suite||!output)die("usage: milpbench --suite FILE --out NEW_DIR [--limit SECONDS --seed N --solvers FILE --checker FILE]");
    Instance *instances=NULL;Solver *solvers=NULL;size_t ni=load_instances(suite,&instances),ns=load_solvers(config,&solvers);if(!ni||!ns)die("empty suite or solver list");
    char *checker=absolute(checker_arg),*self=absolute(argv[0]);char checker_sha[65],runner_sha[65];if(lab_sha256(checker,checker_sha)||lab_sha256(self,runner_sha))die("cannot hash evaluator");
    mkdirs(output);char *out=absolute(output);char *allpath=pathcat(out,"runs.jsonl");FILE *all=fopen(allpath,"wx");if(!all)die("output directory already contains results");
    char *envpath=pathcat(out,"environment.json");FILE *env=fopen(envpath,"wx");if(!env)die("environment output exists");
    struct utsname machine;if(uname(&machine))die("uname");
    fprintf(env,"{\"schema\":1,\"system\":");lab_json_string(env,machine.sysname);fprintf(env,",\"machine\":");lab_json_string(env,machine.machine);fprintf(env,",\"os_release\":");lab_json_string(env,machine.release);
    char h[65];if(lab_sha256(suite,h))die("suite hash");fprintf(env,",\"suite_sha256\":");lab_json_string(env,h);if(lab_sha256(config,h))die("config hash");fprintf(env,",\"solver_config_sha256\":");lab_json_string(env,h);
    fprintf(env,",\"runner_sha256\":");lab_json_string(env,runner_sha);fprintf(env,",\"checker_sha256\":");lab_json_string(env,checker_sha);
    fprintf(env,",\"time_limit_s\":%.9g,\"seed\":%d,\"requested_threads\":1,\"solver_order\":\"rotated_by_instance\",\"timing\":\"end-to-end solver process; checker excluded\",\"argv\":",limit,seed);argv_json(env,argv);fprintf(env,"}\n");fclose(env);
    if(ni>max)ni=max;
    for(size_t i=0;i<ni;i++){verify(instances[i].path,instances[i].sha);for(size_t j=0;j<ns;j++)one(&instances[i],&solvers[(j+i)%ns],checker,checker_sha,runner_sha,out,limit,grace,seed,all);}
    fclose(all);for(size_t i=0;i<ni;i++){free(instances[i].id);free(instances[i].path);free(instances[i].sha);free(instances[i].reference);free(instances[i].value);}free(instances);
    for(size_t j=0;j<ns;j++){free(solvers[j].name);free(solvers[j].binary);free(solvers[j].sha);free(solvers[j].version);}free(solvers);free(checker);free(self);free(out);free(allpath);free(envpath);return 0;
}
