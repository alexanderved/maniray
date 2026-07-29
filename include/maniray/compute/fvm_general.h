#ifndef _MR_FVM_GENERAL_H
#define _MR_FVM_GENERAL_H

#include "maniray/utils/types.h"
#include "maniray/compute/octree.h"

typedef enum mr_node_type {
    MR_NODE_TYPE_NONE,
    MR_NODE_TYPE_BOUNDARY,
    MR_NODE_TYPE_INTERPOLATION,
    MR_NODE_TYPE_OTHER,
} mr_node_type;

typedef struct mr_discretization_data {
    mr_node_type type;

    union {
        mr_int grid_connection; // NONE
        mr_int donor_node_idx; // INTERPOLATION
    };
} mr_discretization_data;

void mr_fit_grids_to_charts(mr_ocforest *forest);
void mr_connect_overset_grids(mr_ocforest *forest);

#endif // _MR_FVM_GENERAL_H