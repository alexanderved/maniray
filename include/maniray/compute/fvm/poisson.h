#ifndef _MR_FVM_POISSON_H
#define _MR_FVM_POISSON_H

#include "maniray/compute/manifold.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/codes.h"
#include "maniray/compute/fvm/boundary.h"
#include "maniray/compute/linear_system/matrix.h"
#include "maniray/compute/linear_system/vector.h"

#define MR_POISSON_SOLUTION_EXTRA_FIELD 1

typedef struct mr_fvm_poisson_solution {
    mr_float value;
} mr_fvm_poisson_solution;

typedef struct mr_fvm_poisson mr_fvm_poisson;
typedef mr_float (*mr_fvm_poisson_source_fn)(mr_fvm_poisson *poisson, mr_int cell_idx);

struct mr_fvm_poisson {
    mr_ocforest *forest;
    mr_code_map *code_map;

    mr_fvm_poisson_source_fn source_fn;
    mr_boundary_condition *bc;

    mr_sparse_matrix *discr_mat;
    mr_vector *source_terms;
};

mr_fvm_poisson *mr_fvm_poisson_create();
void mr_fvm_poisson_destroy(mr_fvm_poisson *poisson);

mr_ocforest *mr_fvm_poisson_ocforest_initialize(
    mr_fvm_poisson *poisson,
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots
);
mr_ocforest *mr_fvm_poisson_ocforest_update(mr_fvm_poisson *poisson);
void mr_fvm_poisson_ocforest_finalize(mr_fvm_poisson *poisson);

void mr_fvm_poisson_set_source_term_fn(mr_fvm_poisson *poisson, mr_fvm_poisson_source_fn source_fn);
void mr_fvm_poisson_set_boundary_condition(mr_fvm_poisson *poisson, mr_boundary_condition *bc);

int mr_fvm_poisson_build_discretization_matrix(mr_fvm_poisson *poisson);
int mr_fvm_poisson_build_source_terms(mr_fvm_poisson *poisson);
int mr_fvm_poisson_solve(mr_fvm_poisson *poisson);

#endif // _MR_FVM_POISSON_H