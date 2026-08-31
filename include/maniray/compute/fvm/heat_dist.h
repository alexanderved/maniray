#ifndef _MR_FVM_HEAT_DIST_H
#define _MR_FVM_HEAT_DIST_H

#include "maniray/compute/manifold.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/codes.h"
#include "maniray/compute/fvm/boundary.h"
#include "maniray/compute/linear_system/solver.h"

#define MR_HEAT_DIST_SOLUTION_EXTRA_FIELD 1

typedef struct mr_fvm_heat_distance_solution {
    mr_float grad[MR_NB_AXES];
    mr_float dist;
} mr_fvm_heat_distance_solution;

typedef struct mr_fvm_heat_distance mr_fvm_heat_distance;
typedef mr_float (*mr_fvm_heat_distance_initial_cond_fn)(mr_fvm_heat_distance *heat_distance, mr_int cell_idx);
typedef mr_sign (*mr_fvm_heat_distance_sign_fn)(mr_fvm_heat_distance *heat_distance, mr_uint chart_idx, mr_float p[MR_NB_AXES]);

struct mr_fvm_heat_distance {
    mr_linear_system_solver *solver;

    mr_ocforest *forest;
    mr_code_map *code_map;

    mr_fvm_heat_distance_initial_cond_fn init_cond_fn;
    mr_fvm_heat_distance_sign_fn sign_fn;
    mr_boundary_condition *bc;

    mr_sparse_matrix *heat_eq_mat;
    mr_sparse_matrix *poisson_mat;

    mr_vector *init_cond_terms;
};

mr_fvm_heat_distance *mr_fvm_heat_distance_create();
void mr_fvm_heat_distance_destroy(mr_fvm_heat_distance *heat_distance);

mr_ocforest *mr_fvm_heat_distance_ocforest_initialize(
    mr_fvm_heat_distance *heat_distance,
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots
);
mr_ocforest *mr_fvm_heat_distance_ocforest_update(mr_fvm_heat_distance *heat_distance);
void mr_fvm_heat_distance_ocforest_finalize(mr_fvm_heat_distance *heat_distance);

void mr_fvm_heat_distance_set_boundary_condition(mr_fvm_heat_distance *heat_distance, mr_boundary_condition *bc);
void mr_fvm_heat_distance_set_initial_condition_fn(mr_fvm_heat_distance *heat_distance, mr_fvm_heat_distance_initial_cond_fn init_cond_fn);
void mr_fvm_heat_distance_set_sign_fn(mr_fvm_heat_distance *heat_distance, mr_fvm_heat_distance_sign_fn sign_fn);

int mr_fvm_heat_distance_build_discretization_matrix(mr_fvm_heat_distance *heat_distance);
int mr_fvm_heat_distance_build_initial_condition_terms(mr_fvm_heat_distance *heat_distance);

int mr_fvm_heat_distance_solve(mr_fvm_heat_distance *heat_distance);

#endif // _MR_FVM_HEAT_DIST_H