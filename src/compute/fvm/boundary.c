#include <assert.h>
#include <tgmath.h>

#include "maniray/compute/fvm/boundary.h"
#include "maniray/compute/fvm/grid.h"

#define GET_COORD(var, elem, axis) \
    do { \
        switch (axis) { \
            case MR_AXIS_X: var = (elem)->x; break; \
            case MR_AXIS_Y: var = (elem)->y; break; \
            case MR_AXIS_Z: var = (elem)->z; break; \
            default: assert(false); \
        } \
    } while (0)

#define EPSILON (1e-6f)

bool mr_is_boundary_cell(mr_ocforest *forest, mr_int cell_idx, mr_direction dir) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type != MR_CELL_TYPE_BOUNDARY) {
        return false;
    }

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);

    mr_int local_idx = cell_idx - node->first_child;
    if (!mr_is_cell_local_idx_face_adjacent(local_idx, dir)) {
        return false;
    }

    mr_axis axis = mr_direction_get_axis(dir);
    mr_float sign_mul = mr_direction_get_sign_mul(dir);

    mr_octree_root *root = mr_ocforest_get_root(forest, node->root);
    if (root->flags & (1 << axis)) {
        return false;
    }
    mr_octree_node *root_node = mr_ocforest_get_node(forest, root->node_idx);

    mr_float cell_coord, root_coord;
    GET_COORD(cell_coord, cell, axis);
    GET_COORD(root_coord, root_node, axis);

    mr_float cell_face = cell_coord + sign_mul * cell->dim / 2.0f;
    mr_float root_face = root_coord + sign_mul * root_node->dim / 2.0f;

    return fabs(cell_face - root_face) <= EPSILON;
}