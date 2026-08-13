#ifndef _MR_FVM_INTERPOLATION_H
#define _MR_FVM_INTERPOLATION_H

#include "maniray/compute/octree.h"

typedef int (*mr_fvm_interpolation_fn)(mr_ocforest *forest, mr_int node_idx, mr_float coef, void *userdata);
typedef struct mr_fvm_interpolation_cb {
    mr_fvm_interpolation_fn fn;
    void *userdata;
} mr_fvm_interpolation_cb;

static inline mr_fvm_interpolation_cb mr_fvm_interpolation_cb_create(mr_fvm_interpolation_fn fn, void *userdata) {
    return (mr_fvm_interpolation_cb) { fn, userdata };
}

static inline mr_fvm_interpolation_cb mr_fvm_interpolation_cb_null() {
    return (mr_fvm_interpolation_cb) { NULL, NULL };
}

bool mr_fvm_point_has_interpolation_stencil(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES]);
int mr_fvm_perform_interpolation(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES], mr_fvm_interpolation_cb interp);

bool mr_fvm_admits_ghost_cell(mr_ocforest *forest, mr_int node_idx, mr_direction dir[MR_ADJACENCY_VERTEX]);
int mr_fvm_calculate_ghost_cell(mr_ocforest *forest, mr_int node_idx, mr_direction dir[MR_ADJACENCY_VERTEX], mr_fvm_interpolation_cb interp);

#endif // _MR_FVM_INTERPOLATION_H