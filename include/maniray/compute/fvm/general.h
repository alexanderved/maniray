#ifndef _MR_FVM_GENERAL_H
#define _MR_FVM_GENERAL_H

#include "maniray/utils/types.h"
#include "maniray/compute/octree.h"

#define MR_DISCR_DATA_EXTRA_FIELD 0

typedef enum mr_cell_type {
    MR_CELL_TYPE_NONE,
    MR_CELL_TYPE_BOUNDARY,
    MR_CELL_TYPE_INTERPOLATION,
    MR_CELL_TYPE_OTHER,
} mr_cell_type;

typedef struct mr_discretization_data {
    mr_cell_type type;

    union {
        mr_int grid_connection; // NONE
        mr_int donor_root_idx;  // INTERPOLATION
    };
} mr_discretization_data;

/*
 * TODO:
 * 1. Add refinement on chart boundaries
 * 2. Add refinement for interface interpolation (might need to implement FLAG_REFINE and refine_flaged function)
*/

void mr_fvm_fit_grids_to_charts(mr_ocforest *forest);
void mr_fvm_connect_overset_grids(mr_ocforest *forest);

#endif // _MR_FVM_GENERAL_H