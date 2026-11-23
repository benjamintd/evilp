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
#include <sys/select.h>
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
    char configuration_sha[65], overlay_sha[65];
} Solver;
typedef struct {
    int code, signal, timed_out;
    double wall, cpu, suspended;
    long rss;
} Process;
static volatile sig_atomic_t interrupted = 0;
static void child_exited(int signo) {
    (void)signo;
}
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
/* A second clock records host suspension independently of solver deadlines.
 * Read it inside the awake-time interval so descheduling cannot fake sleep. */
static double including_suspend(void) {
#if defined(__APPLE__) || defined(CLOCK_BOOTTIME)
    struct timespec ts;
#ifdef __APPLE__
    clockid_t id = CLOCK_MONOTONIC_RAW;
#else
    clockid_t id = CLOCK_BOOTTIME;
#endif
    if (clock_gettime(id, &ts) == 0)
        return (double)ts.tv_sec + (double)ts.tv_nsec / 1e9;
#endif
    return -1;
}
static Process run(char *const args[], const char *log, double limit) {
    Process r = {0};
    double begin = lab_now();
    double including_begin = including_suspend();
    /* Block notifications across waitid and atomically unblock them in pselect.
     * This avoids charging short-lived solvers the old polling interval. */
    sigset_t blocked, previous, waiting;
    sigemptyset(&blocked);
    sigaddset(&blocked, SIGCHLD);
    sigaddset(&blocked, SIGINT);
    sigaddset(&blocked, SIGTERM);
    if (sigprocmask(SIG_BLOCK, &blocked, &previous))
        die("sigprocmask");
    waiting = previous;
    sigdelset(&waiting, SIGCHLD);
    sigdelset(&waiting, SIGINT);
    sigdelset(&waiting, SIGTERM);
    pid_t pid = fork();
    if (pid < 0)
        die("fork");
    if (!pid) {
        if (sigprocmask(SIG_SETMASK, &previous, NULL) || setpgid(0, 0))
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
        double until = terminated ? kill_at : limit;
        double delay = terminated == 2 ? 1 : fmax(0, fmin(1, until - elapsed));
        struct timespec pause = {(time_t)delay, (long)((delay - (time_t)delay) * 1e9)};
        if (pselect(0, NULL, NULL, NULL, &pause, &waiting) < 0 && errno != EINTR)
            die("pselect");
    }
    if (sigprocmask(SIG_SETMASK, &previous, NULL))
        die("sigprocmask restore");
    double including_end = including_suspend();
    r.wall = lab_now() - begin;
    r.suspended = including_begin >= 0 && including_end >= 0
                      ? fmax(0.0, including_end - including_begin - r.wall) : -1;
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
    while (**p == ' ' || **p == '\t' || **p == '\n' || **p == '\r')
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
static const char *const evidence_fields[] = {"search_route", "closure_reason", "bound_origin", "bound_scope"};
enum { EVIDENCE_FIELDS = sizeof evidence_fields / sizeof *evidence_fields };
typedef struct {
    char value[EVIDENCE_FIELDS][64];
} SearchEvidence;

static const char *json_status(const char *p, const char *const known[], size_t count, char *bound, SearchEvidence *evidence) {
    json_space(&p);
    if (*p++ != '{')
        return NULL;
    char value[64] = "", numeric_bound[128] = "";
    SearchEvidence parsed = {0};
    int seen = 0, seen_bound = 0;
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
        } else if (bound && !strcmp(key, "best_bound")) {
            const char *start = p;
            if (seen_bound++ || json_value(&p, 0))
                return NULL;
            if (*start == '-' || isdigit((unsigned char)*start)) {
                size_t n = (size_t)(p - start);
                if (n >= sizeof numeric_bound)
                    return NULL;
                memcpy(numeric_bound, start, n);
                if (!isfinite(strtod(numeric_bound, NULL)))
                    return NULL;
            } else if (p - start != 4 || memcmp(start, "null", 4)) {
                return NULL;
            }
        } else {
            size_t field = 0;
            while (field < EVIDENCE_FIELDS && strcmp(key, evidence_fields[field]))
                ++field;
            if (evidence && field < EVIDENCE_FIELDS) {
                if (parsed.value[field][0] || json_string(&p, parsed.value[field], sizeof parsed.value[field]) ||
                    !parsed.value[field][0])
                    return NULL;
            } else if (json_value(&p, 0)) {
                return NULL;
            }
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
    for (size_t i = 0; i < count; ++i)
        if (!strcmp(value, known[i])) {
            if (bound)
                strcpy(bound, numeric_bound);
            if (evidence)
                *evidence = parsed;
            return known[i];
        }
    return NULL;
}
static const char *lab_status(const char *p, char bound[128], SearchEvidence *evidence) {
    static const char *const known[] = {"optimal",     "infeasible",  "unbounded",  "numerical_failure",
                                        "unsupported", "input_error", "time_limit", "iteration_limit",
                                        "node_limit",  "gap_limit", "lp_unbounded"};
    return json_status(p, known, sizeof known / sizeof *known, bound, evidence);
}

/* milpcheck returns a structured result for each of its five normal exit codes. */
static const char *checker_report_status(const char *report, int code) {
    static const char *const statuses[] = {"valid", "invalid", "error", "unsupported", "no_solution",
                                           "unavailable"};
    const char *status = json_status(report, statuses, sizeof statuses / sizeof *statuses, NULL, NULL);
    if (code < 0 || code > 4 || !status)
        return NULL;
    return !strcmp(status, statuses[code]) || (code == 4 && !strcmp(status, "unavailable"))
               ? status : NULL;
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
/* Cbc sometimes omits a terminal result block after preprocessing, but writes
 * a more specific status in the first line of its native solution file. This
 * is another solver report, never an independent feasibility/proof result. */
static const char *cbc_solution_status(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f)
        return NULL;
    char line[16384];
    int have = fgets(line, sizeof line, f) != NULL;
    int bad = ferror(f);
    fclose(f);
    if (!have || bad || (strlen(line) == sizeof line - 1 && !strchr(line, '\n')))
        return NULL;
    static const char *header[] = {"Optimal", "Infeasible", "Integer infeasible", "Unbounded"};
    static const char *status[] = {"optimal", "infeasible", "infeasible", "unbounded"};
    for (size_t i = 0; i < sizeof header / sizeof *header; i++) {
        size_t n = strlen(header[i]);
        if (strncmp(line, header[i], n))
            continue;
        const char *tail = line + n;
        if (!strncmp(tail, " - objective value ", 19)) {
            tail += 19;
            char *end;
            double objective = strtod(tail, &end);
            if (end == tail || !isfinite(objective))
                return NULL;
            tail = end;
        }
        while (isspace((unsigned char)*tail))
            tail++;
        if (!*tail)
            return status[i];
    }
    return NULL;
}
static int cbc_report_conflict(const char *log_status, const char *solution_status) {
    if (!solution_status)
        return 0;
    if (!strcmp(log_status, "lp_unbounded"))
        return strcmp(solution_status, "unbounded") != 0;
    if (strcmp(log_status, "optimal") && strcmp(log_status, "infeasible") &&
        strcmp(log_status, "unbounded"))
        return 0;
    return strcmp(log_status, solution_status) != 0;
}
static const char *status_from_log(const char *solver, const char *path, Process p,
                                   const char **reported, char bound[128], SearchEvidence *evidence) {
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
            const char *parsed = lab_status(line, bound, evidence);
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
            else if (!strncmp(line, "Pre-processing says infeasible or unbounded",
                              sizeof("Pre-processing says infeasible or unbounded") - 1) &&
                     (line[sizeof("Pre-processing says infeasible or unbounded") - 1] == 0 ||
                      isspace((unsigned char)
                                  line[sizeof("Pre-processing says infeasible or unbounded") - 1])))
                status = "infeasible_or_unbounded";
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
        a[n] = (Solver){lab_strdup(v[0]), absolute(v[1]), lab_strdup(v[2]), lab_strdup(v[3]), {0}};
        verify(a[n].binary, a[n].sha);
        n++;
    }
    free(line);
    fclose(f);
    *out = a;
    return n;
}
/* Configuration identity excludes per-instance output paths. The runner hash
 * binds adapter defaults; each attempt retains its actual argv/settings. */
static void freeze_config(Solver *solver, const char *out, const char *runner_sha,
                          int relax, double limit, double grace, int seed, int paired) {
    char name[64];
    solver->overlay_sha[0] = 0;
    snprintf(name, sizeof name, "%s-config.json", solver->name);
    char *path = pathcat(out, name);
    FILE *f = fopen(path, "wx");
    if (!f)
        die("cannot create frozen solver configuration");
    fprintf(f, "{\"runner_sha256\":\"%s\",\"solver_sha256\":\"%s\",\"adapter\":\"%s\","
            "\"mode\":\"%s\",\"time_limit_s\":%.17g,\"watchdog_grace_s\":%.17g,"
            "\"requested_threads\":1,\"seed\":%d,\"paired_caps\":%s",
            runner_sha, solver->sha, solver->name, relax ? "relax" : "mip", limit, grace,
            seed, paired ? "true" : "false");
    if (!strcmp(solver->name, "lab")) {
        const char *overlay = getenv("LAB_SCIP_SET");
        if (overlay) {
            char *copy = pathcat(out, "lab-scip.set");
            FILE *src = fopen(overlay, "rb");
            FILE *dst = fopen(copy, "wbx");
            if (!src || !dst)
                die("cannot snapshot LAB_SCIP_SET");
            char buffer[8192];
            size_t n;
            while ((n = fread(buffer, 1, sizeof buffer, src)))
                if (fwrite(buffer, 1, n, dst) != n)
                    die("writing LAB_SCIP_SET snapshot");
            if (ferror(src) || fclose(src) || fclose(dst) || lab_sha256(copy, solver->overlay_sha) ||
                setenv("LAB_SCIP_SET", copy, 1))
                die("freezing LAB_SCIP_SET");
            free(copy);
        }
        fprintf(f, ",\"scip_set_sha256\":");
        if (*solver->overlay_sha)
            lab_json_string(f, solver->overlay_sha);
        else
            fputs("null", f);
        fputs(",\"subscip_fraction_env\":", f);
        const char *fraction = getenv("LAB_SUBSCIP_FRAC");
        if (fraction)
            lab_json_string(f, fraction);
        else
            fputs("null", f);
        fprintf(f, ",\"phase_trace\":%s", getenv("LAB_PHASE_TRACE") ? "true" : "false");
    }
    fputs("}\n", f);
    if (fclose(f) || lab_sha256(path, solver->configuration_sha))
        die("hashing frozen configuration");
    free(path);
}
static void paired_fields(FILE *f, const char *campaign_id, const char *host, double limit,
                          double grace, int seed) {
    fputs(",\"campaign_id\":", f);
    lab_json_string(f, campaign_id);
    fputs(",\"hardware_id\":", f);
    lab_json_string(f, host);
    fprintf(f, ",\"measurement\":\"solver_process_wall_excluding_checker\","
            "\"budget\":{\"time_limit_s\":%.9g,\"watchdog_grace_s\":%.9g,"
            "\"seed\":%d,\"requested_threads\":1,\"memory_limit_bytes\":null,"
            "\"node_limit\":10000000,\"iteration_limit\":100000}", limit, grace, seed);
}

static void write_campaign(const char *out, const char *id, const char *host, const Solver *solvers,
                           double limit, double grace, int seed, int quiescent,
                           const char *checker_sha, const char *runner_sha, const char *suite_sha) {
    char *path = pathcat(out, "campaign.json");
    FILE *f = fopen(path, "wx");
    if (!f)
        die("cannot create paired campaign");
    fprintf(f, "{\"schema\":1,\"fair_timing_campaign\":%s,\"quiescent\":%s,"
            "\"serial_solver_processes\":true,\"coverage_exploration\":%s,"
            "\"checker_sha256\":\"%s\",\"runner_sha256\":\"%s\",\"suite_sha256\":\"%s\"",
            quiescent ? "true" : "false", quiescent ? "true" : "false",
            quiescent ? "false" : "true", checker_sha, runner_sha, suite_sha);
    paired_fields(f, id, host, limit, grace, seed);
    fputs(",\"engines\":{", f);
    for (size_t j = 0; j < 2; j++) {
        fprintf(f, "%s\"%s\":{\"solver_sha256\":\"%s\",\"configuration_sha256\":\"%s\"}",
                j ? "," : "", solvers[j].name, solvers[j].sha, solvers[j].configuration_sha);
    }
    fputs("}}\n", f);
    if (fclose(f))
        die("writing paired campaign");
    free(path);
}
static void one(const Instance *inst, const Solver *s, const char *checker, const char *checker_sha,
                const char *runner_sha, const char *out, double limit, double grace, int seed,
                int relax, const char *campaign_id, const char *host, FILE *all) {
    char id[1024];
    snprintf(id, sizeof id, "%s/%s/seed-%d", inst->id, s->name, seed);
    char *dir = pathcat(out, id);
    mkdirs(dir);
    char *log = pathcat(dir, "solver.log"), *sol = pathcat(dir, "solution.sol"),
         *options = pathcat(dir, "settings.txt"), *result = pathcat(dir, "run.json"),
         *checklog = pathcat(dir, "checker.json"), *dual = pathcat(dir, "dual.sol"),
         *farkas = pathcat(dir, "infeasible.farkas");
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
            args[k++] = "--farkas";
            args[k++] = farkas;
        }
        fprintf(opt,
                "relaxation=%s\nmode=%s\ntime_limit=%s\nthreads=1\nalgorithm=deterministic\n"
                "seed_argument=none\nreference_objectives=%s\n",
                relax ? "true" : "false", relax ? "relax" : "mip", secs, relax ? "LP" : "MIP");
    } else if (!strcmp(s->name, "highs")) {
        if (campaign_id)
            fputs("mip_max_nodes = 10000000\nsimplex_iteration_limit = 100000\n", opt);
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
    verify(s->binary, s->sha);
    if (s->overlay_sha[0])
        verify(getenv("LAB_SCIP_SET"), s->overlay_sha);
    char bound[128];
    const char *reported;
    SearchEvidence evidence = {0};
    const char *claimed = status_from_log(s->name, log, p, &reported, bound, &evidence);
    const char *solution_status = !strcmp(s->name, "cbc") ? cbc_solution_status(sol) : NULL;
    int status_conflict = cbc_report_conflict(reported, solution_status);
    const char *status_source = !strcmp(claimed, reported) ? "solver_log" : "process";
    if (!strcmp(s->name, "cbc") && !interrupted && p.code == 0 && !p.signal && !p.timed_out) {
        if (status_conflict) {
            claimed = "unknown";
            status_source = "conflicting_solver_reports";
        } else if (solution_status &&
                   (!strcmp(reported, "unknown") || !strcmp(reported, "infeasible_or_unbounded"))) {
            /* A solution-file unbounded claim alone does not distinguish MILP
             * unboundedness from an unbounded continuous relaxation. */
            claimed =
                !strcmp(solution_status, "unbounded") ? "infeasible_or_unbounded" : solution_status;
            if (strcmp(claimed, reported))
                status_source = "solution_header";
        }
    }
    Process cp = {.code = -1};
    int has_solution = access(sol, R_OK) == 0;
    int has_dual = relax && !strcmp(s->name, "lab") && access(dual, R_OK) == 0;
    int has_farkas = relax && !strcmp(s->name, "lab") && access(farkas, R_OK) == 0;
    int has_evidence = has_solution || has_dual || has_farkas;
    char solution_sha[65], dual_sha[65], farkas_sha[65];
    if ((has_solution && lab_sha256(sol, solution_sha)) || (has_dual && lab_sha256(dual, dual_sha)) ||
        (has_farkas && lab_sha256(farkas, farkas_sha)))
        die("cannot hash retained checker evidence");
    char *check = NULL;
    char *ca[16] = {NULL};
    if (has_evidence && !interrupted) {
        ca[0] = (char *)checker;
        ca[1] = "--model";
        ca[2] = inst->path;
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
        if (has_farkas) {
            ca[c++] = "--farkas";
            ca[c++] = farkas;
        }
        ca[c] = NULL;
        cp = run(ca, checklog, 120);
        verify(checker, checker_sha);
        check = readfile(checklog);
    }
    const char *checked_status = check && !cp.timed_out && !cp.signal
                                     ? checker_report_status(check, cp.code) : NULL;
    if (!checked_status) {
        free(check);
        check = lab_strdup(has_evidence
                               ? "{\"status\":\"error\",\"checked_feasible\":false,"
                                 "\"message\":\"checker did not complete with a matching JSON report\"}"
                               : "{\"status\":\"unavailable\",\"checked_feasible\":false}");
    }
    size_t checklen = strlen(check);
    while (checklen && isspace((unsigned char)check[checklen - 1]))
        check[--checklen] = 0;
    /* Raw newlines in a validated JSON object can occur only between tokens.
     * Keep one line per attempt in runs.jsonl even for a pretty-printed report. */
    for (size_t i = 0; i < checklen; ++i)
        if (check[i] == '\r' || check[i] == '\n')
            check[i] = ' ';
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
    fprintf(f, ",\"configuration_sha256\":\"%s\"", s->configuration_sha);
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
    if (campaign_id)
        paired_fields(f, campaign_id, host, limit, grace, seed);
    if (p.suspended >= 0)
        fprintf(f, ",\"host_suspend_s\":%.9g,\"sleep_affected\":%s", p.suspended,
                p.suspended > 5.0 ? "true" : "false");
    fprintf(f, ",\"claimed_status\":");
    lab_json_string(f, claimed);
    fprintf(f, ",\"reported_solver_status\":");
    lab_json_string(f, reported);
    fprintf(f, ",\"reported_solution_status\":");
    if (solution_status)
        lab_json_string(f, solution_status);
    else
        fputs("null", f);
    fprintf(f, ",\"status_source\":");
    lab_json_string(f, status_source);
    fprintf(f, ",\"status_conflict\":%s", status_conflict ? "true" : "false");
    for (size_t field = 0; field < EVIDENCE_FIELDS; ++field) {
        fprintf(f, ",\"%s\":", evidence_fields[field]);
        if (evidence.value[field][0])
            lab_json_string(f, evidence.value[field]);
        else
            fputs("null", f);
    }
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
    fputs(",\"checker_argv\":", f);
    argv_json(f, ca);
    fputs(",\"evidence_sha256\":{", f);
    if (has_solution)
        fprintf(f, "\"solution.sol\":\"%s\"", solution_sha);
    if (has_dual)
        fprintf(f, "%s\"dual.sol\":\"%s\"", has_solution ? "," : "", dual_sha);
    if (has_farkas)
        fprintf(f, "%s\"infeasible.farkas\":\"%s\"", has_solution || has_dual ? "," : "", farkas_sha);
    fputc('}', f);
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
           checked_status ? checked_status : (has_evidence ? "error" : "no_solution"));
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
    free(farkas);
}
int main(int argc, char **argv) {
    if (argc == 2 && !strcmp(argv[1], "--help")) {
        puts("milpbench --suite FILE --out NEW_DIR [--solvers FILE] [--checker FILE]\n"
             "          [--limit SECONDS] [--grace SECONDS] [--seed POSITIVE_INT]\n"
             "          [--max-instances N] [--relax | --mip] [--paired [--quiescent]]\n"
             "--relax currently supports only lab; suite reference objectives must be LP values.\n"
             "--mip permits lab/native mixed runs and keeps integrality checks enabled.\n"
             "Run from the project root. Results are immutable per output directory.");
        return 0;
    }
    const char *suite = NULL, *config = "config/solvers.tsv", *output = NULL,
               *checker_arg = "build/milpcheck";
    double limit = 2, grace = 15;
    int seed = 1, relax = 0, mip = 0, paired = 0, quiescent = 0;
    size_t max = SIZE_MAX;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--quiescent")) {
            quiescent = 1;
            continue;
        }
        if (!strcmp(argv[i], "--paired")) {
            paired = 1;
            continue;
        }
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
    if (quiescent && !paired)
        die("--quiescent requires --paired");
    if (paired && (!mip || ns != 2 || seed != 1 ||
        (strcmp(solvers[0].name, "lab") && strcmp(solvers[1].name, "lab")) ||
        (strcmp(solvers[0].name, "highs") && strcmp(solvers[1].name, "highs"))))
        die("--paired requires --mip, seed 1, and exactly lab and highs");
    if (paired && grace <= 0)
        die("paired campaigns require positive watchdog grace");
    if (paired && (getenv("LAB_SCIP_SET") || getenv("LAB_PHASE_TRACE")))
        die("paired adapter defaults require LAB_SCIP_SET and LAB_PHASE_TRACE unset");
    size_t loaded_instances = ni;
    struct sigaction action = {0};
    action.sa_handler = on_signal;
    sigemptyset(&action.sa_mask);
    if (sigaction(SIGINT, &action, NULL) || sigaction(SIGTERM, &action, NULL))
        die("sigaction");
    action.sa_handler = child_exited;
    action.sa_flags = SA_NOCLDSTOP;
    if (sigaction(SIGCHLD, &action, NULL))
        die("SIGCHLD handler");
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
    for (size_t j = 0; j < ns; j++)
        freeze_config(&solvers[j], out, runner_sha, relax, limit, grace, seed, paired);
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
    char campaign_id[128];
    if (paired) {
        struct timespec started;
        if (clock_gettime(CLOCK_REALTIME, &started))
            die("cannot identify campaign");
        snprintf(campaign_id, sizeof campaign_id, "%s-%ld-%ld-%ld", h, (long)started.tv_sec,
                 started.tv_nsec, (long)getpid());
        write_campaign(out, campaign_id, machine.nodename, solvers, limit, grace, seed, quiescent,
                       checker_sha, runner_sha, h);
    }
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
                grace, seed, relax, paired ? campaign_id : NULL, machine.nodename, all);
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
