#ifndef _MR_LIN_SYS_SOLVER_H
#define _MR_LIN_SYS_SOLVER_H

#include <lis.h>

#include "maniray/compute/linear_system/matrix.h"
#include "maniray/compute/linear_system/vector.h"

typedef enum mr_linear_system_solver_type {
    MR_SOLVER_BICG,
    MR_SOLVER_CGS,
    MR_SOLVER_BICGSTAB,
    MR_SOLVER_GPBICG,
    MR_SOLVER_TFQMR,
    MR_SOLVER_GMRES,
} mr_linear_system_solver_type;

typedef enum mr_linear_system_percon_type {
    MR_PRECON_NONE,
    MR_PRECON_JACOBI,
    MR_PRECON_ILU,
    MR_PRECON_SSOR,
    MR_PRECON_ILUT,
    MR_PRECON_SAAMG,
} mr_linear_system_percon_type;

typedef struct mr_linear_system_solver {
    LIS_SOLVER inner;
    LIS_VECTOR solution;

    mr_sparse_matrix *matrix;
} mr_linear_system_solver;

mr_linear_system_solver *mr_linear_system_solver_create();
void mr_linear_system_solver_destroy(mr_linear_system_solver *solver);

int mr_linear_system_solver_print_debug_info(mr_linear_system_solver *solver);

int mr_linear_system_solver_set_options(
    mr_linear_system_solver *solver,
    mr_linear_system_solver_type type,
    mr_linear_system_percon_type precon_type,
    mr_float64 tol
);
int mr_linear_system_solver_set_matrix(mr_linear_system_solver *solver, mr_sparse_matrix *matrix);

int mr_linear_system_solve(mr_linear_system_solver *solver, mr_vector *rhs);
int mr_linear_system_solver_get_solution(mr_linear_system_solver *solver, mr_int i, mr_float64 *value);

#endif // _MR_LIN_SYS_SOLVER_H