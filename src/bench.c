#define _DEFAULT_SOURCE
#define _DARWIN_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#include "common.h"
#include <ctype.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <zlib.h>

typedef struct {
    char *id, *path, *sha, *reference, *value;
} Instance;
typedef struct {
    char *name, *binary, *sha, *version;
} Solver;
typedef struct {
    int code, signal, timed_out;
    double wall, cpu;
    long rss;
} Process;
static volatile sig_atomic_t interrupted = 0;
static void on_signal(int signo) {
    interrupted = signo;
}
static void die(const char *s) {
    fprintf(stderr, "milpbench: %s\n", s);
    exit(2);
}
static char *pathcat(const char *a, const char *b) {
    size_t n = strlen(a) + strlen(b) + 2;
    char *p = lab_alloc(n, 1);
    snprintf(p, n, "%s/%s", a, b);
    return p;
}
static char *absolute(const char *s) {
    char *p = realpath(s, NULL);
    if (!p)
        die(s);
    return p;
}
static void mkdirs(const char *path) {
    char *s = lab_strdup(path);
    for (char *p = s + 1; *p; p++)
        if (*p == '/') {
            *p = 0;
            if (mkdir(s, 0755) && errno != EEXIST)
                die("mkdir");
            *p = '/';
        }
    if (mkdir(s, 0755) && errno != EEXIST)
        die("mkdir");
    free(s);
}
static int safe_id(const char *s) {
    if (!*s || !strcmp(s, ".") || !strcmp(s, ".."))
        return 0;
    for (; *s; s++)
        if (!isalnum((unsigned char)*s) && !strchr("._-", *s))
            return 0;
    return 1;
}
static void verify(const char *path, const char *expected) {
    char h[65];
    if (strlen(expected) != 64 || lab_sha256(path, h) || strcmp(h, expected))
        die("SHA-256 mismatch or unreadable locked file");
}
/* Cbc's MPS reader ignores OBJSENSE: pass the original sense explicitly. */
static int objective_max(const char *path) {
    gzFile f = gzopen(path, "rb");
    if (!f)
        die("cannot read objective sense");
    char line[16384], a[128], b[128];
    int waiting = 0, max = 0;
    while (gzgets(f, line, sizeof line)) {
        if (line[0] == '*')
            continue;
        int n = sscanf(line, "%127s %127s", a, b);
        if (n < 1)
            continue;
        if (!strcmp(a, "ROWS"))
            break;
        if (!strcmp(a, "OBJSENSE")) {
            waiting = 1;
            if (n < 2)
                continue;
            snprintf(a, sizeof a, "%s", b);
        }
        if (waiting) {
            if (!strcmp(a, "MAX") || !strcmp(a, "MAXIMIZE"))
                max = 1;
            else if (strcmp(a, "MIN") && strcmp(a, "MINIMIZE"))
                die("unsupported objective sense");
            break;
        }
    }
    gzclose(f);
    return max;
}
static Process run(char *const args[], const char *log, double limit) {
    Process r = {0};
    double begin = lab_now();
    pid_t pid = fork();
    if (pid < 0)
        die("fork");
    if (!pid) {
        if (setpgid(0, 0))
            _exit(126);
        int fd = open(log, O_WRONLY | O_CREAT | O_EXCL, 0644);
        if (fd < 0)
            _exit(126);
        if (dup2(fd, STDOUT_FILENO) < 0 || dup2(fd, STDERR_FILENO) < 0)
            _exit(126);
        close(fd);
        int input = open("/dev/null", O_RDONLY);
        if (input >= 0) {
            dup2(input, STDIN_FILENO);
            close(input);
        }
        setenv("OMP_NUM_THREADS", "1", 1);
        setenv("OPENBLAS_NUM_THREADS", "1", 1);
        setenv("MKL_NUM_THREADS", "1", 1);
        execv(args[0], args);
        perror("execv");
        _exit(127);
    }
    (void)setpgid(pid, pid);
    int status = 0;
    struct rusage usage;
    memset(&usage, 0, sizeof usage);
    int terminated = 0;
    double kill_at = 0;
    for (;;) {
        siginfo_t info = {0};
        if (waitid(P_PID, (id_t)pid, &info, WEXITED | WNOHANG | WNOWAIT) < 0) {
            if (errno == EINTR)
                continue;
            die("waitid");
        }
        if (info.si_pid == pid) {
            /* Keep the PID reserved until group cleanup, then collect resources. */
            kill(-pid, SIGKILL);
            pid_t reaped;
            do {
                reaped = wait4(pid, &status, 0, &usage);
            } while (reaped < 0 && errno == EINTR);
            if (reaped != pid)
                die("wait4");
            break;
        }
        double elapsed = lab_now() - begin;
        if (!terminated && (elapsed > limit || interrupted)) {
            r.timed_out = interrupted ? 0 : 1;
            terminated = 1;
            kill_at = elapsed + 1;
            kill(-pid, SIGTERM);
        }
        if (terminated == 1 && elapsed > kill_at) {
            kill(-pid, SIGKILL);
            terminated = 2;
        }
        struct timespec pause = {0, 20000000};
        nanosleep(&pause, NULL);
    }
    r.wall = lab_now() - begin;
    r.code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
    r.signal = WIFSIGNALED(status) ? WTERMSIG(status) : 0;
    r.cpu = (double)usage.ru_utime.tv_sec + (double)usage.ru_utime.tv_usec / 1e6 +
            (double)usage.ru_stime.tv_sec + (double)usage.ru_stime.tv_usec / 1e6;
#ifdef __APPLE__
    r.rss = usage.ru_maxrss / 1024;
#else
    r.rss = usage.ru_maxrss;
#endif
    return r;
}
static char *readfile(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f)
        return lab_strdup("");
    if (fseek(f, 0, SEEK_END)) {
        fclose(f);
        return lab_strdup("");
    }
    long n = ftell(f);
    if (n < 0 || n > 16 * 1024 * 1024)
        die("oversized result");
    rewind(f);
    char *s = lab_alloc((size_t)n + 1, 1);
    size_t got = fread(s, 1, (size_t)n, f);
    s[got] = 0;
    fclose(f);
    return s;
}
static void argv_json(FILE *f, char *const args[]) {
    fputc('[', f);
    for (int i = 0; args[i]; i++) {
        if (i)
            fputc(',', f);
        lab_json_string(f, args[i]);
    }
    fputc(']', f);
}
/* The lab adapter emits one JSON result object. Parse its top-level status,
 * rather than matching words inside an explanatory message. Other values are
 * structurally validated but do not become independent numerical evidence. */
static void json_space(const char **p) {
    while (isspace((unsigned char)**p))
        ++*p;
}
static int json_string(const char **p, char *out, size_t cap) {
    if (**p != '"')
        return -1;
    ++*p;
    size_t n = 0;
    while (**p && **p != '"') {
        unsigned c = (unsigned char)*(*p)++;
        if (c < 32)
            return -1;
        if (c == '\\') {
            c = (unsigned char)*(*p)++;
            if (!c)
                return -1;
            if (c == 'u') {
                c = 0;
                for (int i = 0; i < 4; ++i) {
                    unsigned h = (unsigned char)**p;
                    if (!isxdigit((unsigned char)h))
                        return -1;
                    c = 16 * c + (h <= '9' ? h - '0' : (unsigned)tolower((int)h) - 'a' + 10);
                    ++*p;
                }
            } else if (strchr("\"\\/", (int)c)) {
                /* Literal escaped punctuation. */
            } else if (strchr("bfnrt", (int)c)) {
                c = 128;
            } else {
                return -1;
            }
        }
        /* Status names and keys are ASCII. Preserve non-ASCII/control escapes
         * as a nonmatching byte, never an embedded C-string terminator. */
        if (out && n + 1 < cap)
            out[n] = (char)(c >= 32 && c < 127 ? c : 128);
        ++n;
    }
    if (**p != '"')
        return -1;
    ++*p;
    if (out && cap) {
        if (n >= cap)
            out[0] = 0;
        else
            out[n] = 0;
    }
    return 0;
}
static int json_value(const char **p, unsigned depth) {
    if (depth > 64)
        return -1;
    json_space(p);
    if (**p == '"')
        return json_string(p, NULL, 0);
    if (**p == '{' || **p == '[') {
        int object = **p == '{';
        char end = object ? '}' : ']';
        ++*p;
        json_space(p);
        if (**p == end) {
            ++*p;
            return 0;
        }
        for (;;) {
            if (object) {
                if (json_string(p, NULL, 0))
                    return -1;
                json_space(p);
                if (*(*p)++ != ':')
                    return -1;
            }
            if (json_value(p, depth + 1))
                return -1;
            json_space(p);
            if (**p == end) {
                ++*p;
                return 0;
            }
            if (*(*p)++ != ',')
                return -1;
            json_space(p);
        }
    }
    static const char *literals[] = {"true", "false", "null"};
    for (size_t i = 0; i < 3; ++i) {
        size_t n = strlen(literals[i]);
        if (!strncmp(*p, literals[i], n)) {
            *p += n;
            return 0;
        }
    }
    if (**p == '-')
        ++*p;
    if (**p == '0')
        ++*p;
    else {
        if (**p < '1' || **p > '9')
            return -1;
        while (isdigit((unsigned char)**p))
            ++*p;
    }
    if (**p == '.') {
        ++*p;
        if (!isdigit((unsigned char)**p))
            return -1;
        while (isdigit((unsigned char)**p))
            ++*p;
    }
    if (**p == 'e' || **p == 'E') {
        ++*p;
        if (**p == '-' || **p == '+')
            ++*p;
        if (!isdigit((unsigned char)**p))
            return -1;
        while (isdigit((unsigned char)**p))
            ++*p;
    }
    return 0;
}
static const char *lab_status(const char *p) {
    static const char *known[] = {"optimal",     "infeasible",  "unbounded",  "numerical_failure",
                                  "unsupported", "input_error", "time_limit", "iteration_limit",
                                  "node_limit",  "lp_unbounded"};
    json_space(&p);
    if (*p++ != '{')
        return NULL;
    char value[64] = "";
    int seen = 0;
    json_space(&p);
    if (*p == '}')
        return NULL;
    for (;;) {
        char key[64];
        if (json_string(&p, key, sizeof key))
            return NULL;
        json_space(&p);
        if (*p++ != ':')
            return NULL;
        json_space(&p);
        if (!strcmp(key, "status")) {
            if (seen++ || json_string(&p, value, sizeof value))
                return NULL;
        } else if (json_value(&p, 0)) {
            return NULL;
        }
        json_space(&p);
        if (*p == '}') {
            ++p;
            break;
        }
        if (*p++ != ',')
            return NULL;
        json_space(&p);
    }
    json_space(&p);
    if (*p || !seen)
        return NULL;
    for (size_t i = 0; i < sizeof known / sizeof *known; ++i)
        if (!strcmp(value, known[i]))
            return known[i];
    return NULL;
}
static int cbc_lp_optimal_line(const char *line) {
    const char *value;
    int summary = 0;
    if (!strncmp(line, "Optimal - objective value ", 26))
        value = line + 26;
    else if (!strncmp(line, "Optimal objective ", 18)) {
        value = line + 18;
        summary = 1;
    } else {
        return 0;
    }
    char *end;
    double objective = strtod(value, &end);
    if (end == value || !isfinite(objective))
        return 0;
    if (summary)
        return !strncmp(end, " - ", 3) && strstr(end, " iterations time ") != NULL;
    while (isspace((unsigned char)*end))
        ++end;
    return !*end;
}
static const char *status_from_log(const char *solver, const char *path, Process p,
                                   const char **reported, char bound[128]) {
    *reported = "unknown";
    bound[0] = 0;
    FILE *f = fopen(path, "r");
    if (!f)
        return "launch_error";
    char line[16384];
    const char *status = "unknown";
    int highs_report = 0, input_timeout = 0, cbc_lp_optimal = 0, cbc_mip = 0;
    while (fgets(line, sizeof line, f)) {
        if (!strcmp(solver, "lab")) {
            size_t length = strlen(line);
            if (length == sizeof line - 1 && line[length - 1] != '\n') {
                int c;
                while ((c = fgetc(f)) != EOF && c != '\n') {
                }
                continue;
            }
            const char *parsed = lab_status(line);
            if (parsed)
                status = parsed;
            continue;
        }
        if (!strcmp(solver, "highs") && strstr(line, "Solving report"))
            highs_report = 1;
        if (!strcmp(solver, "highs") && strstr(line, "Parser reached timeout"))
            input_timeout = 1;
        if (!strcmp(solver, "highs") &&
            (strstr(line, "Model status") || (highs_report && !strncmp(line, "  Status ", 9)))) {
            if (strstr(line, "Infeasible or unbounded") ||
                strstr(line, "Unbounded or infeasible") ||
                strstr(line, "infeasible or unbounded") || strstr(line, "unbounded or infeasible"))
                status = "infeasible_or_unbounded";
            else if (strstr(line, "Optimal"))
                status = "optimal";
            else if (strstr(line, "Infeasible"))
                status = "infeasible";
            else if (strstr(line, "Unbounded"))
                status = "unbounded";
            else if (strstr(line, "Time limit"))
                status = "time_limit";
        } else if (!strcmp(solver, "scip") && strstr(line, "SCIP Status")) {
            if (strstr(line, "infeasible or unbounded") || strstr(line, "unbounded or infeasible"))
                status = "infeasible_or_unbounded";
            else if (strstr(line, "optimal solution"))
                status = "optimal";
            else if (strstr(line, "infeasible"))
                status = "infeasible";
            else if (strstr(line, "unbounded"))
                status = "unbounded";
            else if (strstr(line, "time limit"))
                status = "time_limit";
        } else if (!strcmp(solver, "cbc")) {
            if ((!strncmp(line, "Cbc", 3) && isdigit((unsigned char)line[3])) ||
                !strncmp(line, "Continuous objective value is ", 30))
                cbc_mip = 1;
            if (cbc_lp_optimal_line(line))
                cbc_lp_optimal = 1;
            if (strstr(line, "Result - Optimal solution found"))
                status = "optimal";
            else if (strstr(line, "Result - Problem proven infeasible") ||
                     strstr(line, "Problem is infeasible"))
                status = "infeasible";
            else if (strstr(line, "Result - Linear relaxation unbounded"))
                status = "lp_unbounded";
            else if (strstr(line, "Result - Stopped on time limit"))
                status = "time_limit";
        }
        char *start = NULL;
        if (!strcmp(solver, "highs"))
            start = strstr(line, "Dual bound");
        else if (!strcmp(solver, "scip"))
            start = strstr(line, "Dual Bound");
        else
            start = strstr(line, "Lower bound:");
        if (start) {
            while (*start && !isdigit((unsigned char)*start) && *start != '-' && *start != '+')
                start++;
            char *end;
            double d = strtod(start, &end);
            if (end != start && isfinite(d))
                snprintf(bound, 128, "%.17g", d);
        }
    }
    fclose(f);
    /* A continuous-only Cbc solve uses Clp's terminal output, without Cbc's
     * MILP result block. Do not promote a root LP message from a MIP run. */
    if (!strcmp(solver, "cbc") && !strcmp(status, "unknown") && cbc_lp_optimal && !cbc_mip)
        status = "optimal";
    /* Preserve the solver's report independently of interruption, watchdog or
     * exit-code classification. A successful-looking report cannot erase a
     * failed process outcome; neither is an independent proof. */
    *reported = status;
    if (interrupted)
        return "interrupted";
    if (p.timed_out)
        return "watchdog_timeout";
    if (p.code == 127 || p.code == 126)
        return "launch_error";
    if (input_timeout)
        return "input_time_limit";
    if (p.signal)
        return "process_error";
    if (!strcmp(solver, "lab") && strcmp(status, "unknown"))
        return status;
    if (p.code != 0)
        return "process_error";
    return status;
}
static int split(char *s, char **v, int max) {
    size_t n = strlen(s);
    while (n && (s[n - 1] == '\n' || s[n - 1] == '\r'))
        s[--n] = 0;
    int count = 0;
    char *p = s;
    for (;;) {
        if (count == max)
            return -1;
        v[count++] = p;
        char *t = strchr(p, '\t');
        if (!t)
            break;
        *t = 0;
        p = t + 1;
    }
    return count;
}
static size_t load_instances(const char *path, Instance **out) {
    FILE *f = fopen(path, "r");
    if (!f)
        die("cannot open suite");
    char *line = NULL, *v[5];
    size_t size = 0, n = 0;
    ssize_t len;
    Instance *a = NULL;
    while ((len = getline(&line, &size, f)) >= 0) {
        if (len <= 1 || line[0] == '#' || !strncmp(line, "id\t", 3))
            continue;
        if (split(line, v, 5) != 5 || !safe_id(v[0]))
            die("malformed suite");
        Instance *next = realloc(a, (n + 1) * sizeof *a);
        if (!next)
            die("allocation");
        a = next;
        for (size_t j = 0; j < n; j++)
            if (!strcmp(a[j].id, v[0]))
                die("duplicate instance ID");
        a[n] = (Instance){lab_strdup(v[0]), absolute(v[1]), lab_strdup(v[2]), lab_strdup(v[3]),
                          lab_strdup(v[4])};
        n++;
    }
    free(line);
    fclose(f);
    *out = a;
    return n;
}
static size_t load_solvers(const char *path, Solver **out) {
    FILE *f = fopen(path, "r");
    if (!f)
        die("cannot open solvers config");
    char *line = NULL, *v[4];
    size_t size = 0, n = 0;
    Solver *a = NULL;
    while (getline(&line, &size, f) >= 0) {
        if (line[0] == '#' || !strncmp(line, "solver\t", 7) || line[0] == '\n')
            continue;
        if (split(line, v, 4) != 4)
            die("malformed solvers config");
        if (strcmp(v[0], "highs") && strcmp(v[0], "scip") && strcmp(v[0], "cbc") &&
            strcmp(v[0], "lab"))
            die("unsupported solver adapter");
        Solver *next = realloc(a, (n + 1) * sizeof *a);
        if (!next)
            die("allocation");
        a = next;
        for (size_t j = 0; j < n; j++)
            if (!strcmp(a[j].name, v[0]))
                die("duplicate solver");
        a[n] = (Solver){lab_strdup(v[0]), absolute(v[1]), lab_strdup(v[2]), lab_strdup(v[3])};
        verify(a[n].binary, a[n].sha);
        n++;
    }
    free(line);
    fclose(f);
    *out = a;
    return n;
}
static void one(const Instance *inst, const Solver *s, const char *checker, const char *checker_sha,
                const char *runner_sha, const char *out, double limit, double grace, int seed,
                int relax, FILE *all) {
    char id[1024];
    snprintf(id, sizeof id, "%s/%s/seed-%d", inst->id, s->name, seed);
    char *dir = pathcat(out, id);
    mkdirs(dir);
    char *log = pathcat(dir, "solver.log"), *sol = pathcat(dir, "solution.sol"),
         *options = pathcat(dir, "settings.txt"), *result = pathcat(dir, "run.json"),
         *checklog = pathcat(dir, "checker.json"), *dual = pathcat(dir, "dual.sol");
    char secs[64], seedstr[32];
    snprintf(secs, sizeof secs, "%.9g", limit);
    snprintf(seedstr, sizeof seedstr, "%d", seed);
    char *args[40] = {0};
    int k = 0;
    args[k++] = s->binary;
    char commands[11][PATH_MAX + 128];
    memset(commands, 0, sizeof commands);
    FILE *opt = fopen(options, "wx");
    if (!opt)
        die("result directory is not new");
    if (!strcmp(s->name, "lab")) {
        args[k++] = inst->path;
        args[k++] = relax ? "--relax" : "--mip";
        args[k++] = "--time-limit";
        args[k++] = secs;
        args[k++] = "--solution";
        args[k++] = sol;
        if (relax) {
            args[k++] = "--dual";
            args[k++] = dual;
        }
        fprintf(opt,
                "relaxation=%s\nmode=%s\ntime_limit=%s\nthreads=1\nalgorithm=deterministic\n"
                "seed_argument=none\nreference_objectives=%s\n",
                relax ? "true" : "false", relax ? "relax" : "mip", secs, relax ? "LP" : "MIP");
    } else if (!strcmp(s->name, "highs")) {
        fprintf(opt,
                "time_limit = %s\nthreads = 1\nrandom_seed = %s\nmip_rel_gap = 0\nmip_abs_gap = "
                "0\nwrite_solution_to_file = true\nwrite_solution_style = 0\nsolution_file = "
                "%s\nlog_file = %s/highs.log\n",
                secs, seedstr, sol, dir);
        args[k++] = inst->path;
        args[k++] = "--options_file";
        args[k++] = options;
    } else if (!strcmp(s->name, "cbc")) {
        args[k++] = inst->path;
        args[k++] = objective_max(inst->path) ? "-max" : "-min";
        args[k++] = "-threads";
        args[k++] = "1";
        args[k++] = "-randomSeed";
        args[k++] = seedstr;
        args[k++] = "-randomCbcSeed";
        args[k++] = seedstr;
        args[k++] = "-timeMode";
        args[k++] = "elapsed";
        args[k++] = "-seconds";
        args[k++] = secs;
        args[k++] = "-ratioGap";
        args[k++] = "0";
        args[k++] = "-allowableGap";
        args[k++] = "0";
        args[k++] = "-solve";
        args[k++] = "-solution";
        args[k++] = sol;
        args[k++] = "-quit";
        fprintf(opt,
                "objective_sense=%s\nthreads=1\nrandomSeed=%s\nrandomCbcSeed=%s\ntimeMode="
                "elapsed\nseconds=%s\nratioGap=0\nallowableGap=0\n",
                args[2], seedstr, seedstr, secs);
    } else {
        if (strpbrk(inst->path, "\"\n\r") || strpbrk(sol, "\"\n\r"))
            die("SCIP adapter cannot represent quotes/newlines in filenames");
        snprintf(commands[0], sizeof commands[0], "set limits time %s", secs);
        snprintf(commands[1], sizeof commands[1], "set parallel maxnthreads 1");
        snprintf(commands[2], sizeof commands[2], "set randomization randomseedshift %s", seedstr);
        snprintf(commands[3], sizeof commands[3], "set limits gap 0");
        snprintf(commands[4], sizeof commands[4], "set limits absgap 0");
        snprintf(commands[5], sizeof commands[5], "read \"%s\"", inst->path);
        snprintf(commands[6], sizeof commands[6], "optimize");
        snprintf(commands[7], sizeof commands[7], "write solution \"%s\"", sol);
        snprintf(commands[8], sizeof commands[8], "quit");
        snprintf(commands[9], sizeof commands[9], "set timing clocktype 2");
        snprintf(commands[10], sizeof commands[10], "set presolving milp threads 1");
        args[k++] = "-c";
        args[k++] = commands[9];
        fprintf(opt, "%s\n", commands[9]);
        args[k++] = "-c";
        args[k++] = commands[10];
        fprintf(opt, "%s\n", commands[10]);
        for (int i = 0; i < 9; i++) {
            args[k++] = "-c";
            args[k++] = commands[i];
            fprintf(opt, "%s\n", commands[i]);
        }
    }
    fclose(opt);
    args[k] = NULL;
    Process p = run(args, log, limit + grace);
    char bound[128];
    const char *reported;
    const char *claimed = status_from_log(s->name, log, p, &reported, bound);
    Process cp = {.code = -1};
    int has_solution = access(sol, R_OK) == 0;
    int has_dual = relax && !strcmp(s->name, "lab") && access(dual, R_OK) == 0;
    int has_evidence = has_solution || has_dual;
    char *check = NULL;
    if (has_evidence && !interrupted) {
        char *ca[16] = {(char *)checker, "--model", inst->path};
        int c = 3;
        if (has_solution) {
            ca[c++] = "--solution";
            ca[c++] = sol;
            ca[c++] = "--format";
            ca[c++] = !strcmp(s->name, "lab") ? "simple" : s->name;
        }
        if (relax)
            ca[c++] = "--relax";
        if (has_dual) {
            ca[c++] = "--dual";
            ca[c++] = dual;
        }
        ca[c] = NULL;
        cp = run(ca, checklog, 120);
        check = readfile(checklog);
    }
    if (!check || check[0] != '{' || cp.timed_out || cp.signal || cp.code < 0 || cp.code > 4) {
        free(check);
        check = lab_strdup("{\"status\":\"unavailable\",\"checked_feasible\":false}");
    }
    size_t checklen = strlen(check);
    while (checklen && isspace((unsigned char)check[checklen - 1]))
        check[--checklen] = 0;
    FILE *f = fopen(result, "wx");
    if (!f)
        die("cannot write run result");
    fprintf(f, "{\"schema\":1,\"instance\":");
    lab_json_string(f, inst->id);
    fprintf(f, ",\"model_path\":");
    lab_json_string(f, inst->path);
    fprintf(f, ",\"model_sha256\":");
    lab_json_string(f, inst->sha);
    fprintf(f, ",\"solver\":");
    lab_json_string(f, s->name);
    fprintf(f, ",\"solver_version\":");
    lab_json_string(f, s->version);
    fprintf(f, ",\"solver_sha256\":");
    lab_json_string(f, s->sha);
    fprintf(f, ",\"checker_sha256\":");
    lab_json_string(f, checker_sha);
    fprintf(f, ",\"runner_sha256\":");
    lab_json_string(f, runner_sha);
    fprintf(f, ",\"relaxation\":%s,\"mode\":\"%s\",\"dual_file\":", relax ? "true" : "false",
            relax ? "relax" : "mip");
    if (has_dual)
        lab_json_string(f, dual);
    else
        fputs("null", f);
    fprintf(f,
            ",\"seed\":%d,\"requested_threads\":1,\"time_limit_s\":%.9g,\"watchdog_grace_s\":%.9g,"
            "\"wall_s\":%.9g,\"cpu_s\":%.9g,\"peak_rss_kib\":%ld,\"exit_code\":%d,\"signal\":%d,"
            "\"watchdog_timeout\":%s",
            seed, limit, grace, p.wall, p.cpu, p.rss, p.code, p.signal,
            p.timed_out ? "true" : "false");
    fprintf(f, ",\"claimed_status\":");
    lab_json_string(f, claimed);
    fprintf(f, ",\"reported_solver_status\":");
    lab_json_string(f, reported);
    fprintf(f, ",\"reported_dual_bound\":");
    if (bound[0])
        lab_json_string(f, bound);
    else
        fputs("null", f);
    fprintf(f, ",\"reference_status\":");
    lab_json_string(f, inst->reference);
    fprintf(f, ",\"reference_objective\":");
    lab_json_string(f, inst->value);
    fprintf(f, ",\"argv\":");
    argv_json(f, args);
    fprintf(f, ",\"checker_exit_code\":");
    if (has_evidence)
        fprintf(f, "%d", cp.code);
    else
        fputs("null", f);
    fprintf(f, ",\"checker_wall_s\":%.9g,\"checker\":%s}\n", cp.wall, check);
    fclose(f);
    char *json = readfile(result);
    fputs(json, all);
    fflush(all);
    printf("%s %s seed=%d status=%s wall=%.3fs check=%s\n", inst->id, s->name, seed, claimed,
           p.wall,
           has_evidence ? (cp.code == 0   ? "valid"
                           : cp.code == 1 ? "invalid"
                                          : "unavailable")
                        : "no_solution");
    fflush(stdout);
    free(json);
    free(check);
    free(dir);
    free(log);
    free(sol);
    free(options);
    free(result);
    free(checklog);
    free(dual);
}
int main(int argc, char **argv) {
    if (argc == 2 && !strcmp(argv[1], "--help")) {
        puts("milpbench --suite FILE --out NEW_DIR [--solvers FILE] [--checker FILE]\n"
             "          [--limit SECONDS] [--grace SECONDS] [--seed POSITIVE_INT]\n"
             "          [--max-instances N] [--relax | --mip]\n"
             "--relax currently supports only lab; suite reference objectives must be LP values.\n"
             "--mip permits lab/native mixed runs and keeps integrality checks enabled.\n"
             "Run from the project root. Results are immutable per output directory.");
        return 0;
    }
    const char *suite = NULL, *config = "config/solvers.tsv", *output = NULL,
               *checker_arg = "build/milpcheck";
    double limit = 2, grace = 15;
    int seed = 1, relax = 0, mip = 0;
    size_t max = SIZE_MAX;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--relax")) {
            relax = 1;
            continue;
        }
        if (!strcmp(argv[i], "--mip")) {
            mip = 1;
            continue;
        }
        if (i + 1 >= argc)
            die("missing argument value");
        const char *key = argv[i], *value = argv[++i];
        char *end;
        if (!strcmp(key, "--suite"))
            suite = value;
        else if (!strcmp(key, "--solvers"))
            config = value;
        else if (!strcmp(key, "--out"))
            output = value;
        else if (!strcmp(key, "--checker"))
            checker_arg = value;
        else if (!strcmp(key, "--limit")) {
            limit = strtod(value, &end);
            if (*end || !isfinite(limit) || limit <= 0 || limit > 86400)
                die("invalid time limit");
        } else if (!strcmp(key, "--grace")) {
            grace = strtod(value, &end);
            if (*end || !isfinite(grace) || grace < 0 || grace > 3600)
                die("invalid grace");
        } else if (!strcmp(key, "--seed")) {
            long x = strtol(value, &end, 10);
            if (*end || x < 1 || x > INT_MAX)
                die("seed must be 1..INT_MAX; zero is wall-clock seeded in some solvers");
            seed = (int)x;
        } else if (!strcmp(key, "--max-instances")) {
            unsigned long x = strtoul(value, &end, 10);
            if (*end || !x || value[0] == '-')
                die("invalid instance count");
            max = (size_t)x;
        } else
            die("unknown option");
    }
    if (!suite || !output)
        die("usage: milpbench --suite FILE --out NEW_DIR [--limit SECONDS --seed N --solvers FILE "
            "--checker FILE]");
    if (relax && mip)
        die("--relax and --mip are mutually exclusive");
    Instance *instances = NULL;
    Solver *solvers = NULL;
    size_t ni = load_instances(suite, &instances), ns = load_solvers(config, &solvers);
    if (!ni || !ns)
        die("empty suite or solver list");
    for (size_t j = 0; j < ns; ++j) {
        if (!strcmp(solvers[j].name, "lab") && !relax && !mip)
            die("lab adapter requires --relax or --mip; suite reference objectives must match the "
                "mode");
        if (relax && strcmp(solvers[j].name, "lab"))
            die("--relax currently supports only lab; native MILP adapters cannot be mixed with LP "
                "references");
    }
    size_t loaded_instances = ni;
    struct sigaction action = {0};
    action.sa_handler = on_signal;
    sigemptyset(&action.sa_mask);
    if (sigaction(SIGINT, &action, NULL) || sigaction(SIGTERM, &action, NULL))
        die("sigaction");
    char *checker = absolute(checker_arg), *self = absolute(argv[0]);
    char checker_sha[65], runner_sha[65];
    if (lab_sha256(checker, checker_sha) || lab_sha256(self, runner_sha))
        die("cannot hash evaluator");
    mkdirs(output);
    char *out = absolute(output);
    char *allpath = pathcat(out, "runs.jsonl");
    FILE *all = fopen(allpath, "wx");
    if (!all)
        die("output directory already contains results");
    char *envpath = pathcat(out, "environment.json");
    FILE *env = fopen(envpath, "wx");
    if (!env)
        die("environment output exists");
    struct utsname machine;
    if (uname(&machine))
        die("uname");
    fprintf(env, "{\"schema\":1,\"system\":");
    lab_json_string(env, machine.sysname);
    fprintf(env, ",\"machine\":");
    lab_json_string(env, machine.machine);
    fprintf(env, ",\"os_release\":");
    lab_json_string(env, machine.release);
    char h[65];
    if (lab_sha256(suite, h))
        die("suite hash");
    fprintf(env, ",\"suite_sha256\":");
    lab_json_string(env, h);
    if (lab_sha256(config, h))
        die("config hash");
    fprintf(env, ",\"solver_config_sha256\":");
    lab_json_string(env, h);
    fprintf(env, ",\"runner_sha256\":");
    lab_json_string(env, runner_sha);
    fprintf(env, ",\"checker_sha256\":");
    lab_json_string(env, checker_sha);
    fprintf(env,
            ",\"relaxation\":%s,\"mode\":\"%s\",\"time_limit_s\":%.9g,\"seed\":%d,\"requested_"
            "threads\":1,\"solver_"
            "order\":\"rotated_"
            "by_instance\",\"timing\":\"end-to-end solver process; checker excluded\",\"argv\":",
            relax ? "true" : "false", relax ? "relax" : "mip", limit, seed);
    argv_json(env, argv);
    fprintf(env, "}\n");
    fclose(env);
    if (ni > max)
        ni = max;
    for (size_t i = 0; i < ni && !interrupted; i++) {
        verify(instances[i].path, instances[i].sha);
        for (size_t j = 0; j < ns && !interrupted; j++)
            one(&instances[i], &solvers[(j + i) % ns], checker, checker_sha, runner_sha, out, limit,
                grace, seed, relax, all);
    }
    fclose(all);
    for (size_t i = 0; i < loaded_instances; i++) {
        free(instances[i].id);
        free(instances[i].path);
        free(instances[i].sha);
        free(instances[i].reference);
        free(instances[i].value);
    }
    free(instances);
    for (size_t j = 0; j < ns; j++) {
        free(solvers[j].name);
        free(solvers[j].binary);
        free(solvers[j].sha);
        free(solvers[j].version);
    }
    free(solvers);
    free(checker);
    free(self);
    free(out);
    free(allpath);
    free(envpath);
    return interrupted ? 128 + interrupted : 0;
}
