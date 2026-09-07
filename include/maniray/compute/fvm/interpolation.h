#ifndef _MR_FVM_INTERPOLATION_H
#define _MR_FVM_INTERPOLATION_H

#include "maniray/compute/octree.h"

MR_DEFINE_CALLBACK(mr_fvm_interpolation, int, mr_ocforest *forest, mr_int cell_idx, mr_float coef)

int mr_fvm_perform_interpolation(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES], mr_fvm_interpolation_cb interp);
int mr_fvm_calculate_ghost_cell(mr_ocforest *forest, mr_int coarse_cell_idx, mr_int fine_cell_idx, mr_fvm_interpolation_cb interp);
int mr_fvm_interpolate_face_value(mr_ocforest *forest, mr_int cell_idx, mr_direction face_dir, mr_fvm_interpolation_cb interp);

/*
 * TODO:
 * 1. Add choice for order of accuracy
*/

#endif // _MR_FVM_INTERPOLATION_H