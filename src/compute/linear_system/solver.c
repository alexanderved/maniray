#include <assert.h>

#include "maniray/compute/linear_system/solver.h"
#include "maniray/utils/xmalloc.h"
#include "maniray/utils/misc.h"

mr_linear_system_solver *mr_linear_system_solver_create() {
    mr_linear_system_solver *solver = xmalloc(sizeof(mr_linear_system_solver));

    if (lis_solver_create(&solver->inner) != LIS_SUCCESS) {
        free(solver);

        return NULL;
    }

    solver->solution = NULL;
    solver->matrix = NULL;

    return solver;
}

void mr_linear_system_solver_destroy(mr_linear_system_solver *solver) {
    if (!solver) {
        return;
    }

    lis_vector_destroy(solver->solution);
    lis_solver_destroy(solver->inner);
}

int mr_linear_system_solver_print_debug_info(mr_linear_system_solver *solver) {
    return lis_solver_set_option("-print 2", solver->inner) == LIS_SUCCESS ? MR_SUCCESS : MR_FAILURE;
}

// TODO: Make list of types larger and add more customizations
int mr_linear_system_solver_set_options(
    mr_linear_system_solver *solver,
    mr_linear_system_solver_type type,
    mr_linear_system_percon_type precon_type,
    mr_float64 tol
) {
    assert(solver);

    static const LIS_INT type_map[] = {
        [MR_SOLVER_BICG] = LIS_SOLVER_BICG,
        [MR_SOLVER_CGS] = LIS_SOLVER_CGS,
        [MR_SOLVER_BICGSTAB] = LIS_SOLVER_BICGSTAB,
        [MR_SOLVER_GPBICG] = LIS_SOLVER_GPBICG,
        [MR_SOLVER_TFQMR] = LIS_SOLVER_TFQMR,
        [MR_SOLVER_GMRES] = LIS_SOLVER_GMRES,
    };

    static const LIS_INT precon_type_map[] = {
        [MR_PRECON_NONE] = LIS_PRECON_TYPE_NONE,
        [MR_PRECON_JACOBI] = LIS_PRECON_TYPE_JACOBI,
        [MR_PRECON_ILU] = LIS_PRECON_TYPE_ILU,
        [MR_PRECON_SSOR] = LIS_PRECON_TYPE_SSOR,
        [MR_PRECON_ILUT] = LIS_PRECON_TYPE_ILUT,
        [MR_PRECON_SAAMG] = LIS_PRECON_TYPE_SAAMG,
    };

    solver->inner->options[LIS_OPTIONS_SOLVER] = type_map[type];
    solver->inner->options[LIS_OPTIONS_PRECON] = precon_type_map[precon_type];
    solver->inner->params[LIS_PARAMS_RESID - LIS_OPTIONS_LEN] = tol;

    return MR_SUCCESS;
}

int mr_linear_system_solver_set_matrix(mr_linear_system_solver *solver, mr_sparse_matrix *matrix) {
    assert(solver);
    assert(matrix);

    if (solver->solution) {
        lis_vector_destroy(solver->solution);
        solver->solution = NULL;
    }

    if (lis_vector_duplicate(matrix->inner, &solver->solution) != LIS_SUCCESS) {
        return MR_FAILURE;
    }

    solver->matrix = matrix;

    return MR_SUCCESS;
}

// TODO: Reuse preconditioners
int mr_linear_system_solve(mr_linear_system_solver *solver, mr_vector *rhs) {
    assert(solver);
    assert(rhs);
    assert(solver->matrix);

    return lis_solve(solver->matrix->inner, rhs->inner, solver->solution, solver->inner) == LIS_SUCCESS ? MR_SUCCESS
                                                                                                        : MR_FAILURE;
}

int mr_linear_system_solver_get_solution(mr_linear_system_solver *solver, mr_int i, mr_float64 *value) {
    assert(solver);
    assert(value);
    assert(solver->solution);

    return lis_vector_get_value(solver->solution, i, value) == LIS_SUCCESS ? MR_SUCCESS : MR_FAILURE;
}