#ifndef LAB_SEARCH_RESULT_H
#define LAB_SEARCH_RESULT_H

/* A search observation in the objective units of the model passed to that
 * search. A returned point borrows the caller's buffer. Closure is a solver
 * claim; only the independent checker can certify evidence. */
typedef enum {
    LAB_MIP_CONTINUE, LAB_MIP_OPTIMAL, LAB_MIP_INFEASIBLE, LAB_MIP_TIME_LIMIT,
    LAB_MIP_NODE_LIMIT, LAB_MIP_ITERATION_LIMIT, LAB_MIP_GAP_LIMIT,
    LAB_MIP_NUMERICAL_FAILURE, LAB_MIP_UNSUPPORTED, LAB_MIP_LP_UNBOUNDED, LAB_MIP_INPUT_ERROR
} LabMipStatus;

typedef enum {
    LAB_CLOSE_NONE, LAB_CLOSE_BACKEND, LAB_CLOSE_EXHAUSTED, LAB_CLOSE_MATCHING_BOUND,
    LAB_CLOSE_INTEGER_LATTICE, LAB_CLOSE_UNSAT, LAB_CLOSE_DOMAIN, LAB_CLOSE_UNIQUE_POINT
} LabClosure;

typedef enum {
    LAB_BOUND_NONE, LAB_BOUND_BACKEND, LAB_BOUND_EXHAUSTION, LAB_BOUND_STRUCTURAL,
    LAB_BOUND_INTEGER_LATTICE, LAB_BOUND_FRONTIER, LAB_BOUND_UNIQUE_POINT
} LabBoundOrigin;

/* The scope of a present search bound is the whole MIP supplied to the route.
 * Node-local LP bounds stay in the tree. A caller lifting a result to another
 * model must establish equivalence in original objective units before publishing. */
typedef struct {
    LabBoundOrigin origin;
    double value;
} LabSearchBound;

typedef struct {
    const char *route;
    LabMipStatus status;
    LabClosure closure;
    const double *point;
    LabSearchBound bound;
} LabSearchResult;

const char *lab_mip_status_name(LabMipStatus status);
const char *lab_closure_name(LabClosure closure);
const char *lab_bound_origin_name(LabBoundOrigin origin);
#endif
