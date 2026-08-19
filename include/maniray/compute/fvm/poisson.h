#ifndef _MR_FVM_POISSON_H
#define _MR_FVM_POISSON_H

#include "maniray/compute/manifold.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/codes.h"

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
};

mr_fvm_poisson *mr_fvm_poisson_create(
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots,
    mr_fvm_poisson_source_fn source_fn
);
void mr_fvm_poisson_destroy(mr_fvm_poisson *poisson);

mr_ocforest *mr_fvm_poisson_get_ocforest(mr_fvm_poisson *poisson);

int mr_fvm_poisson_build_discretization_matrix(mr_fvm_poisson *poisson);
int mr_fvm_poisson_build_source_terms(mr_fvm_poisson *poisson)

#endif // _MR_FVM_POISSON_H