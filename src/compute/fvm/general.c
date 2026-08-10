#include <stdio.h>
#include <string.h>

#include "maniray/compute/fvm/general.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/utils/misc.h"

static bool is_node_in_bounds(mr_ocforest *forest, size_t chart_idx, mr_octree_node *node) {
    mr_float center[MR_NB_AXES] = { node->x, node->y, node->z };

    if (!mr_manifold_is_in_bounds(forest->manifold, chart_idx, center)) {
        return false;
    }

    for (size_t i = 0; i < 8; ++i) {
        mr_float hdim = node->dim / 2.0f;
        mr_float x = center[0] + hdim * (i & 1 ? 1.0f : -1.0f);
        mr_float y = center[1] + hdim * ((i >> 1) & 1 ? 1.0f : -1.0f);
        mr_float z = center[2] + hdim * ((i >> 2) & 1 ? 1.0f : -1.0f);

        mr_float p[MR_NB_AXES] = { x, y, z };
        if (!mr_manifold_is_in_bounds(forest->manifold, chart_idx, p)) {
            return false;
        }
    }

    return true;
}

static int activate_in_bounds(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (is_node_in_bounds(forest, node->chart_idx, node)) {
        node->flags |= MR_OCTREE_NODE_FLAG_ACTIVE;
    } else {
        node->flags &= ~MR_OCTREE_NODE_FLAG_ACTIVE;
        discr_data->type = MR_CELL_TYPE_OTHER;
    }

    return MR_SUCCESS;
}

void mr_fvm_fit_grids_to_charts(mr_ocforest *forest) {
    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_leaves_apply(forest, octree_idx, mr_octree_cond_cb_null(), mr_octree_apply_cb_create(activate_in_bounds, NULL), false);
    }
}

static int setup_initial_conds(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type == MR_CELL_TYPE_NONE) {
        discr_data->grid_connection = forest->nb_roots;
    }

    return MR_SUCCESS;
}

static int transition_center(mr_ocforest *forest, mr_int idx, mr_int other_octree_idx, mr_float res[MR_NB_AXES]) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    mr_int other_root_idx = forest->roots[other_octree_idx].node_idx;
    mr_octree_node *other_root_node = mr_ocforest_get_node(forest, other_root_idx);

    mr_float orig_coords[] = { node->x, node->y, node->z };
    mr_float new_coords[] = { 0.0f, 0.0f, 0.0f };
    if (mr_manifold_transition(forest->manifold, node->chart_idx, other_root_node->chart_idx, new_coords, orig_coords) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    memcpy(res, new_coords, MR_NB_AXES * sizeof(mr_float));

    return MR_SUCCESS;
}

static mr_int locate_center_in_other_octree(mr_ocforest *forest, mr_int idx, mr_int other_octree_idx) {
    mr_float coords[] = { 0.0f, 0.0f, 0.0f };
    if (transition_center(forest, idx, other_octree_idx, coords) != MR_SUCCESS) {
        return MR_INVALID_INDEX;
    }

    return mr_octree_locate_point(forest, other_octree_idx, coords);
}

static int disable_nodes_near_boundaries(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type != MR_CELL_TYPE_BOUNDARY) {
        return MR_SUCCESS;
    }

    for (mr_index other_octree_idx = 0; (size_t)other_octree_idx < forest->nb_roots; ++other_octree_idx) {
        if (other_octree_idx == node->root) {
            continue;
        }

        mr_int disabled_node_idx = locate_center_in_other_octree(forest, idx, other_octree_idx);
        if (disabled_node_idx == MR_INVALID_INDEX) {
            continue;
        }

        mr_octree_node *disabled_node = mr_ocforest_get_node(forest, disabled_node_idx);
        mr_discretization_data *disabled_discr_data = mr_ocforest_get_extra(forest, disabled_node_idx, MR_DISCR_DATA_EXTRA_FIELD);

        disabled_node->flags &= ~MR_OCTREE_NODE_FLAG_ACTIVE;
        disabled_discr_data->type = MR_CELL_TYPE_OTHER;
    }

    return MR_SUCCESS;
}

static bool is_valid_interior_node(mr_ocforest *forest, mr_int idx) {
    for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
        mr_int nidx = mr_octree_find_face_neighbor(forest, idx, dir);
        mr_octree_node *n = mr_ocforest_get_node(forest, nidx);
        if (!n) {
            return false;
        }

        if (n->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            if (!(n->flags & MR_OCTREE_NODE_FLAG_ACTIVE)) {
                return false;
            } else {
                continue;
            }
        }

        for (mr_int local_idx = 0; local_idx < MR_OCTREE_NB_CHILDREN; ++local_idx) {
            if ((local_idx >> (dir / 2) & 1) == (dir % 2)) {
                continue;
            }

            mr_octree_node *child = mr_ocforest_get_node(forest, n->first_child + local_idx);
            if (!(child->flags & MR_OCTREE_NODE_FLAG_ACTIVE)) {
                return false;
            }
        }
    }

    return true;
}

static bool is_valid_interpolation_node(mr_ocforest *forest, mr_int idx, mr_int other_octree_idx) {
    mr_float coords[] = { 0.0f, 0.0f, 0.0f };
    if (transition_center(forest, idx, other_octree_idx, coords) != MR_SUCCESS) {
        return false;
    }

    return mr_fvm_point_has_interpolation_stencil(forest, other_octree_idx, coords);
}

static int setup_interpolation(mr_ocforest *forest, mr_int idx, void *userdata) {
    size_t *nb_updates = userdata;

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    mr_int conn = discr_data->grid_connection;
    bool is_valid = false;
    do {
        if (conn == node->root + 1) {
            if (is_valid_interior_node(forest, idx)) {
                is_valid = true;
            } else {
                --discr_data->grid_connection;
                ++(*nb_updates);
            }
        } else if (is_valid_interpolation_node(forest, idx, conn - 1)) {
            is_valid = true;
        } else {
            --discr_data->grid_connection;
            ++(*nb_updates);
        }

        --conn;
    } while (conn != 0 && !is_valid);

    if (discr_data->grid_connection == 0) {
        node->flags &= ~MR_OCTREE_NODE_FLAG_ACTIVE;
        discr_data->type = MR_CELL_TYPE_OTHER;
    }

    return MR_SUCCESS;
}

void mr_fvm_connect_overset_grids(mr_ocforest *forest) {
    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_leaves_apply(forest, octree_idx, mr_octree_cond_cb_null(), mr_octree_apply_cb_create(setup_initial_conds, NULL), false);
    }

    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_leaves_apply(
            forest,
            octree_idx,
            mr_octree_cond_cb_null(),
            mr_octree_apply_cb_create(disable_nodes_near_boundaries, NULL),
            false
        );
    }

    size_t nb_updates = 0;
    do {
        nb_updates = 0;

        for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
            mr_octree_leaves_apply(
                forest,
                octree_idx,
                mr_octree_cond_cb_null(),
                mr_octree_apply_cb_create(setup_interpolation, &nb_updates),
                false
            );
        }
    } while (nb_updates != 0);
}