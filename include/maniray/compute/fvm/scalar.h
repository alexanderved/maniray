#ifndef _MR_FVM_SCALAR_H
#define _MR_FVM_SCALAR_H

#include "maniray/compute/octree.h"
#include "maniray/compute/fvm/boundary.h"

MR_DEFINE_CALLBACK(mr_fvm_scalar_store_coef, int, mr_float coef)

int mr_fvm_scalar_interpolate(mr_ocforest *forest, mr_int idx, mr_fvm_scalar_store_coef_cb cb);

int mr_fvm_scalar_calc_boundary_flux_implicit_term(
    mr_ocforest *forest,
    mr_int idx,
    mr_boundary_condition *cond,
    mr_direction dir,
    mr_fvm_scalar_store_coef_cb cb
);
int mr_fvm_scalar_calc_boundary_flux_source_term(
    mr_ocforest *forest,
    mr_int idx,
    mr_boundary_condition *cond,
    mr_direction dir,
    mr_fvm_scalar_store_coef_cb cb
);

int mr_fvm_scalar_calc_internal_flux(mr_ocforest *forest, mr_int idx, mr_direction dir, mr_fvm_scalar_store_coef_cb cb);

#endif // _MR_FVM_SCALAR_H