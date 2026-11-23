#include "search_result.h"

const char *lab_mip_status_name(LabMipStatus status) {
    static const char *const names[] = {
        "continue", "optimal", "infeasible", "time_limit", "node_limit", "iteration_limit",
        "gap_limit", "numerical_failure", "unsupported", "lp_unbounded", "input_error"
    };
    return names[status];
}

const char *lab_closure_name(LabClosure closure) {
    static const char *const names[] = {
        "none", "backend_status", "search_exhausted", "matching_bound",
        "integer_lattice", "unsatisfiable", "empty_domain", "unique_point"
    };
    return names[closure];
}

const char *lab_bound_origin_name(LabBoundOrigin origin) {
    static const char *const names[] = {
        "none", "backend_dual", "search_exhaustion", "structural_bound",
        "integer_lattice", "open_nodes", "unique_point",
        "objective_box", "source_equalities", "row_penalty", "sat_cutoff_unsat"
    };
    return names[origin];
}
