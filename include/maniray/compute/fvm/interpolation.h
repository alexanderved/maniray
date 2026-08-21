#ifndef _MR_FVM_INTERPOLATION_H
#define _MR_FVM_INTERPOLATION_H

#include "maniray/compute/octree.h"

MR_DEFINE_CALLBACK(mr_fvm_interpolation, int, mr_ocforest *forest, mr_int node_idx, mr_float coef)

int mr_fvm_perform_interpolation(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES], mr_fvm_interpolation_cb interp);
int mr_fvm_calculate_ghost_cell(mr_ocforest *forest, mr_int node_idx, mr_direction dir[MR_ADJACENCY_VERTEX], mr_fvm_interpolation_cb interp);

#endif // _MR_FVM_INTERPOLATION_H