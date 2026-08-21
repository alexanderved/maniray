#include <stdio.h>
#include <string.h>
#include <time.h>

#include "maniray/compute/math.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/utils/misc.h"

static bool is_node_in_bounds(mr_ocforest *forest, mr_int node_idx) {
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    mr_float center[MR_NB_AXES] = { node->x, node->y, node->z };

    if (!mr_manifold_is_in_bounds(forest->manifold, node->chart_idx, center)) {
        return false;
    }

    for (size_t i = 0; i < 8; ++i) {
        mr_float hdim = node->dim / 2.0f;
        mr_float x = center[0] + hdim * (i & 1 ? 1.0f : -1.0f);
        mr_float y = center[1] + hdim * ((i >> 1) & 1 ? 1.0f : -1.0f);
        mr_float z = center[2] + hdim * ((i >> 2) & 1 ? 1.0f : -1.0f);

        mr_float p[MR_NB_AXES] = { x, y, z };
        if (!mr_manifold_is_in_bounds(forest->manifold, node->chart_idx, p)) {
            return false;
        }
    }

    return true;
}

static void disable_node(mr_ocforest *forest, mr_octree_node *node) {
    node->flags &= ~MR_OCTREE_NODE_FLAG_ACTIVE;
    for (mr_int i = 0; i < MR_OCTREE_NB_CELLS_IN_BLOCK; ++i) {
        mr_int cell_idx = node->first_child + i;
        mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

        discr_data->type = MR_CELL_TYPE_EXTERIOR;
        discr_data->grid_connection = 0;
    }
}

static int activate_in_bounds(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    if (is_node_in_bounds(forest, idx)) {
        node->flags |= MR_OCTREE_NODE_FLAG_ACTIVE;
    } else {
        disable_node(forest, node);
    }

    return MR_SUCCESS;
}

void mr_fvm_fit_grids_to_charts(mr_ocforest *forest) {
    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_leaves_apply(forest, octree_idx, mr_octree_apply_cb_create(activate_in_bounds, NULL));
    }
}

static int setup_initial_conds(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type == MR_CELL_TYPE_NONE) {
        discr_data->grid_connection = forest->nb_roots;
    }

    return MR_SUCCESS;
}

static int transition_center(mr_ocforest *forest, mr_int cell_idx, mr_int other_octree_idx, mr_float res[MR_NB_AXES]) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

    mr_int other_root_idx = forest->roots[other_octree_idx].node_idx;
    mr_octree_node *other_root_node = mr_ocforest_get_node(forest, other_root_idx);

    mr_float orig_coords[] = { cell->x, cell->y, cell->z };
    mr_float new_coords[] = { 0.0f, 0.0f, 0.0f };
    if (mr_manifold_transition(forest->manifold, cell->chart_idx, other_root_node->chart_idx, new_coords, orig_coords) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    memcpy(res, new_coords, MR_NB_AXES * sizeof(mr_float));

    return MR_SUCCESS;
}

static mr_int locate_center_in_other_octree_node(mr_ocforest *forest, mr_int cell_idx, mr_int other_octree_idx) {
    mr_float coords[] = { 0.0f, 0.0f, 0.0f };
    if (transition_center(forest, cell_idx, other_octree_idx, coords) != MR_SUCCESS) {
        return MR_INVALID_INDEX;
    }

    return mr_octree_locate_point_in_leaf(forest, other_octree_idx, coords);
}

static int disable_nodes_near_boundaries(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type != MR_CELL_TYPE_BOUNDARY) {
        return MR_SUCCESS;
    }

    for (mr_index other_octree_idx = 0; (size_t)other_octree_idx < forest->nb_roots; ++other_octree_idx) {
        if (other_octree_idx == node->root) {
            continue;
        }

        mr_int disabled_node_idx = locate_center_in_other_octree_node(forest, cell_idx, other_octree_idx);
        if (disabled_node_idx == MR_INVALID_INDEX) {
            continue;
        }

        mr_octree_node *disabled_node = mr_ocforest_get_node(forest, disabled_node_idx);
        disable_node(forest, disabled_node);
    }

    return MR_SUCCESS;
}

static bool is_valid_interior_cell(mr_ocforest *forest, mr_int cell_idx) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);

    mr_int local_idx = cell_idx - node->first_child;

    for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
        if (!mr_is_cell_local_idx_face_adjacent(local_idx, dir)) {
            continue;
        }

        mr_int nidx = mr_octree_find_face_neighbor_node(forest, cell->parent, dir);
        if (nidx == MR_INVALID_INDEX) {
            return false;
        }

        mr_octree_node *n = mr_ocforest_get_node(forest, nidx);
        if (n->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            if (!(n->flags & MR_OCTREE_NODE_FLAG_ACTIVE)) {
                return false;
            }

            continue;
        }

        for (mr_int local_idx = 0; local_idx < MR_OCTREE_NB_CHILDREN; ++local_idx) {
            if (mr_is_node_local_idx_face_adjacent(local_idx, dir)) {
                continue;
            }

            mr_int child_idx = n->first_child + local_idx;
            mr_octree_node *child = mr_ocforest_get_node(forest, child_idx);
            if (!(child->flags & MR_OCTREE_NODE_FLAG_ACTIVE)) {
                return false;
            }
        }
    }

    return true;
}

static bool is_valid_interpolation_node(mr_ocforest *forest, mr_int cell_idx, mr_int other_octree_idx) {
    mr_int node_idx = locate_center_in_other_octree_node(forest, cell_idx, other_octree_idx);
    if (node_idx == MR_INVALID_INDEX) {
        return false;
    }

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    return node->flags & MR_OCTREE_NODE_FLAG_ACTIVE;
}

static int setup_interpolation(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    size_t *nb_updates = userdata;

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    mr_int conn = discr_data->grid_connection;
    do {
        if (conn == node->root + 1) {
            if (is_valid_interior_cell(forest, cell_idx)) {
                break;
            } else {
                --discr_data->grid_connection;
                ++(*nb_updates);
            }
        } else if (is_valid_interpolation_node(forest, cell_idx, conn - 1)) {
            break;
        } else {
            --discr_data->grid_connection;
            ++(*nb_updates);
        }

        --conn;
    } while (conn != 0);

    if (discr_data->grid_connection == 0) {
        disable_node(forest, node);
    }

    return MR_SUCCESS;
}

static void mark_required_cells_in_node(mr_ocforest *forest, mr_int node_idx) {
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    for (mr_int i = 0; i < MR_OCTREE_NB_CELLS_IN_BLOCK; ++i) {
        mr_int cell_idx = node->first_child + i;
        mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
        if (discr_data->type != MR_CELL_TYPE_NONE) {
            continue;
        }

        discr_data->grid_connection = -MR_ABS(discr_data->grid_connection);
    }

    node->value = 1.0f;
}

static int mark_required_nodes(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    mr_int other_octree_idx = MR_ABS(discr_data->grid_connection) - 1;
    if (other_octree_idx < node->root && other_octree_idx >= 0) {
        mr_int interp_node = locate_center_in_other_octree_node(forest, cell_idx, other_octree_idx);
        mark_required_cells_in_node(forest, interp_node);
    }

    return MR_SUCCESS;
}

static bool is_node_unused(mr_ocforest *forest, mr_octree_node *node) {
    bool has_internal_cells = false;
    bool has_unused_interp_cells = false;

    for (mr_int i = 0; i < MR_OCTREE_NB_CELLS_IN_BLOCK; ++i) {
        mr_int cell_idx = node->first_child + i;
        mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
        if (discr_data->type != MR_CELL_TYPE_NONE) {
            continue;
        }

        has_internal_cells = has_internal_cells || discr_data->grid_connection != node->root + 1;
        has_unused_interp_cells = has_unused_interp_cells || (discr_data->grid_connection > 0 && discr_data->grid_connection != node->root + 1);
    }

    return !has_internal_cells && has_unused_interp_cells;
}

static int remove_unused_interp_nodes(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    if (is_node_unused(forest, node)) {
        disable_node(forest, node);
    }

    return MR_SUCCESS;
}

static int make_interp_nodes_internal(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    if (MR_ABS(discr_data->grid_connection) != node->root + 1 && is_valid_interior_cell(forest, cell_idx)) {
        discr_data->grid_connection = node->root + 1;
    }

    return MR_SUCCESS;
}

static int mark_required_higher_nodes(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    mr_int other_octree_idx = MR_ABS(discr_data->grid_connection) - 1;
    if (other_octree_idx > node->root && other_octree_idx >= 0) {
        mr_int interp_node = locate_center_in_other_octree_node(forest, cell_idx, other_octree_idx);
        mark_required_cells_in_node(forest, interp_node);
    }

    return MR_SUCCESS;
}

static int finalize_grids(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    if (MR_ABS(discr_data->grid_connection) == node->root + 1) {
        discr_data->type = MR_CELL_TYPE_INTERIOR;
    } else {
        discr_data->type = MR_CELL_TYPE_INTERPOLATION;
        discr_data->donor_root_idx = MR_ABS(discr_data->grid_connection) - 1;
    }

    return MR_SUCCESS;
}

void mr_fvm_connect_overset_grids(mr_ocforest *forest) {
    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(forest, octree_idx, mr_octree_apply_cb_create(setup_initial_conds, NULL));
    }

    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(
            forest,
            octree_idx,
            mr_octree_apply_cb_create(disable_nodes_near_boundaries, NULL)
        );
    }

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    size_t cnt = 0;
    size_t nb_updates = 0;
    do {
        nb_updates = 0;

        for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
            mr_octree_cells_apply(
                forest,
                octree_idx,
                mr_octree_apply_cb_create(setup_interpolation, &nb_updates)
            );
        }
        ++cnt;
    } while (nb_updates != 0);

    clock_gettime(CLOCK_MONOTONIC, &end);
    long long elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                           (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Setup interpolation (%lu): %.2f ms\n", cnt, (double)elapsed_us / 1000.0);

    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(
            forest,
            octree_idx,
            mr_octree_apply_cb_create(mark_required_nodes, NULL)
        );
    }

    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_leaves_apply(
            forest,
            octree_idx,
            mr_octree_apply_cb_create(remove_unused_interp_nodes, NULL)
        );

        mr_octree_cells_apply(
            forest,
            octree_idx,
            mr_octree_apply_cb_create(make_interp_nodes_internal, NULL)
        );

        mr_octree_cells_apply(
            forest,
            octree_idx,
            mr_octree_apply_cb_create(mark_required_higher_nodes, NULL)
        );
    }

    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(
            forest,
            octree_idx,
            mr_octree_apply_cb_create(finalize_grids, NULL)
        );
    }
}