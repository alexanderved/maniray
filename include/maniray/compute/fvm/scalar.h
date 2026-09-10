#ifndef _MR_FVM_SCALAR_H
#define _MR_FVM_SCALAR_H

#include "maniray/compute/octree.h"
#include "maniray/compute/fvm/boundary.h"

MR_DEFINE_CALLBACK(mr_fvm_scalar_store_coef, int, mr_ocforest *forest, mr_int cell_idx, mr_float64 coef)

int mr_fvm_scalar_calc_center_derivative(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_axis axis,
    mr_fvm_scalar_store_coef_cb store
);

int mr_fvm_scalar_calc_face_derivative(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_direction face_dir,
    mr_axis axis,
    mr_fvm_scalar_store_coef_cb store
);

int mr_fvm_scalar_calc_center_gradient(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_fvm_scalar_store_coef_cb component_store[MR_NB_AXES]
);

int mr_fvm_scalar_calc_face_gradient(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_direction face_dir,
    mr_fvm_scalar_store_coef_cb component_store[MR_NB_AXES]
);

int mr_fvm_scalar_mark_inactive_cell(mr_ocforest *forest, mr_int cell_idx, mr_fvm_scalar_store_coef_cb store);
int mr_fvm_scalar_interpolate(mr_ocforest *forest, mr_int cell_idx, mr_fvm_scalar_store_coef_cb store);

int mr_fvm_scalar_calc_boundary_flux(
    mr_ocforest *forest,
    mr_boundary_condition *bc,
    mr_int cell_idx,
    mr_direction dir,
    mr_fvm_scalar_store_coef_cb store_implicit,
    mr_fvm_scalar_store_coef_cb store_rhs
);

// TODO: Rename to `mr_fvm_scalar_calc_internal_diffusion`, add `store_cross` and rename `store` -> `store_self`
int mr_fvm_scalar_calc_internal_flux(mr_ocforest *forest, mr_int cell_idx, mr_direction dir, mr_fvm_scalar_store_coef_cb store);

int mr_fvm_scalar_calc_transient_term(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_float64 step,
    mr_fvm_scalar_store_coef_cb store_implicit,
    mr_fvm_scalar_store_coef_cb store_rhs
);

#endif // _MR_FVM_SCALAR_H