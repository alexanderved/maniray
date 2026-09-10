#include <assert.h>
#include <tgmath.h>

#include "maniray/compute/geometry.h"
#include "maniray/compute/math.h"
#include "maniray/compute/manifold3d.h"

void mr_cell_face_center(mr_ocforest *forest, mr_int idx, mr_direction face_dir, mr_float face_center[MR_NB_AXES]) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);
    assert(face_center);

    mr_axis axis = mr_direction_get_axis(face_dir);
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, idx);

    face_center[MR_AXIS_X] = cell->x;
    face_center[MR_AXIS_Y] = cell->y;
    face_center[MR_AXIS_Z] = cell->z;

    face_center[axis] += mr_direction_get_sign_mul(face_dir) * cell->dim / 2.0f;
}

mr_float mr_cell_center_metric(mr_ocforest *forest, mr_int idx, size_t i, size_t j) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, idx);
    mr_float center[MR_NB_AXES] = { cell->x, cell->y, cell->z };

    return mr_manifold_metric(forest->manifold, cell->chart_idx, center, i, j);
}

mr_float mr_cell_face_center_metric(mr_ocforest *forest, mr_int idx, mr_direction face_dir, size_t i, size_t j) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, idx);

    mr_float face_center[MR_NB_AXES] = { 0.0f };
    mr_cell_face_center(forest, idx, face_dir, face_center);

    return mr_manifold_metric(forest->manifold, cell->chart_idx, face_center, i, j);
}

mr_float mr_cell_center_inv_metric(mr_ocforest *forest, mr_int idx, size_t i, size_t j) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, idx);
    mr_float center[MR_NB_AXES] = { cell->x, cell->y, cell->z };

    return mr_manifold_inv_metric(forest->manifold, cell->chart_idx, center, i, j);
}

mr_float mr_cell_face_center_inv_metric(mr_ocforest *forest, mr_int idx, mr_direction face_dir, size_t i, size_t j) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, idx);

    mr_float face_center[MR_NB_AXES] = { 0.0f };
    mr_cell_face_center(forest, idx, face_dir, face_center);

    return mr_manifold_inv_metric(forest->manifold, cell->chart_idx, face_center, i, j);
}

mr_float mr_cell_volume(mr_ocforest *forest, mr_int idx) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, idx);
    mr_float center[MR_NB_AXES] = { cell->x, cell->y, cell->z };

    return pow(cell->dim, 3.0f) * sqrt(mr_manifold3d_metric_determinant(forest->manifold, cell->chart_idx, center));
}

mr_float mr_cell_face_area(mr_ocforest *forest, mr_int idx, mr_direction face_dir) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);

    mr_axis axis = mr_direction_get_axis(face_dir);
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, idx);

    mr_float face_center[MR_NB_AXES] = { 0.0f };
    mr_cell_face_center(forest, idx, face_dir, face_center);

    return pow(cell->dim, 2.0f) * sqrt(mr_manifold3d_metric_minor(forest->manifold, cell->chart_idx, face_center, axis, axis));
}