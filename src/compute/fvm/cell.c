#include <tgmath.h>

#include "maniray/compute/fvm/cell.h"
#include "maniray/compute/math.h"
#include "maniray/compute/manifold3d.h"

mr_float mr_cell_volume(mr_ocforest *forest, mr_int idx) {
    if (!forest || idx == MR_INVALID_INDEX) {
        return 0.0f;
    }

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_float center[MR_NB_AXES] = { node->x, node->y, node->z };

    return pow(node->dim, 3.0f) * sqrt(mr_manifold3d_metric_determinant(forest->manifold, node->chart_idx, center));
}

mr_float mr_cell_face_area(mr_ocforest *forest, mr_int idx, mr_direction face_dir) {
    if (!forest || idx == MR_INVALID_INDEX) {
        return 0.0f;
    }

    mr_axis axis = mr_direction_get_axis(face_dir);
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    mr_float hdim = node->dim / 2.0f;
    mr_float face_center[MR_NB_AXES] = { node->x, node->y, node->z };
    face_center[axis] += mr_direction_get_sign_mul(face_dir) * hdim;

    return pow(node->dim, 2.0f) * sqrt(mr_manifold3d_metric_minor(forest->manifold, node->chart_idx, face_center, axis, axis));
}

mr_float mr_cell_neighbor_distance(mr_ocforest *forest, mr_int idx, mr_direction dir) {
    if (!forest || idx == MR_INVALID_INDEX) {
        return 0.0f;
    }

    mr_axis axis = mr_direction_get_axis(dir);
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    mr_float hdim = node->dim / 2.0f;
    mr_float middle[MR_NB_AXES] = { node->x, node->y, node->z };
    middle[axis] += mr_direction_get_sign_mul(dir) * hdim;

    return node->dim * sqrt(mr_manifold_metric(forest->manifold, node->chart_idx, middle, axis, axis));
}