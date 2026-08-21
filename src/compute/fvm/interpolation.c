#include <stdio.h>
#include <string.h>

#include "maniray/compute/math.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/utils/misc.h"

// Q for triQuadratic interpolation
#define Q_STENCIL_DIM 3
#define Q_STENCIL_DIM_SQR (Q_STENCIL_DIM * Q_STENCIL_DIM)

static mr_float lagrange_coef(mr_float x, mr_float stencil[Q_STENCIL_DIM], mr_int idx) {
    mr_float res = 1.0f;
    for (mr_int i = 0; i < Q_STENCIL_DIM; ++i) {
        if (i != idx) {
            res *= (x - stencil[i]) / (stencil[idx] - stencil[i]);
        }
    }

    return res;
}

static mr_float lagrange_coef_with_idx(mr_float x, mr_int first_idx, mr_int idx) {
    mr_float stencil[] = {
        (mr_float)first_idx,
        (mr_float)(first_idx + 1),
        (mr_float)(first_idx + 2),
    };

    return lagrange_coef(x, stencil, idx - first_idx);
}

static mr_float calculate_coef(
    const mr_float p[MR_NB_AXES],
    mr_octree_cell *cell,
    mr_int first_local_idx[MR_NB_AXES],
    mr_int local_idx[MR_NB_AXES]
) {
    if (!p || !cell || !first_local_idx) {
        return 0.0f;
    }

    mr_float local_coords[] = {
        (p[MR_AXIS_X] - cell->x) / cell->dim,
        (p[MR_AXIS_Y] - cell->y) / cell->dim,
        (p[MR_AXIS_Z] - cell->z) / cell->dim,
    };

    mr_float coef = 1.0f;
    for (mr_int i = 0; i < Q_STENCIL_DIM; ++i) {
        coef *= lagrange_coef_with_idx(local_coords[i], first_local_idx[i], local_idx[i]);
    }

    return coef;
}

static mr_int get_first_cell_local_coord(mr_int local_coord) {
    return local_coord / 2;
}

int mr_fvm_perform_interpolation(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES], mr_fvm_interpolation_cb interp) {
    if (!forest || !p || !interp.fn) {
        return MR_FAILURE;
    }

    mr_int host_cell_idx = mr_octree_locate_point_in_cell(forest, octree_idx, p);
    if (host_cell_idx == MR_INVALID_INDEX) {
        return MR_FAILURE;
    }

    mr_octree_cell *host_cell = mr_ocforest_get_cell(forest, host_cell_idx);
    mr_octree_node *host_node = mr_ocforest_get_node(forest, host_cell->parent);
    if (!(host_node->flags & MR_OCTREE_NODE_FLAG_ACTIVE)) {
        return MR_FAILURE;
    }

    mr_int host_local_idx = host_cell_idx - host_node->first_child;
    mr_int host_local_coords[] = {
        host_local_idx & 0x3,
        host_local_idx >> 2 & 0x3,
        host_local_idx >> 4 & 0x3,
    };

    mr_int first_local_coords[] = {
        get_first_cell_local_coord(host_local_coords[0]),
        get_first_cell_local_coord(host_local_coords[1]),
        get_first_cell_local_coord(host_local_coords[2]),
    };

    for (mr_int i = 0; i < Q_STENCIL_DIM; ++i) {
        for (mr_int j = 0; j < Q_STENCIL_DIM; ++j) {
            for (mr_int k = 0; k < Q_STENCIL_DIM; ++k) {
                mr_int local_coords[] = {
                    first_local_coords[0] + i,
                    first_local_coords[1] + j,
                    first_local_coords[2] + k,
                };

                mr_int local_idx = (local_coords[2] << 4) + (local_coords[1] << 2) + local_coords[0];
                mr_int cell_idx = host_node->first_child + local_idx;
                mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

                mr_float coef = calculate_coef(p, cell, first_local_coords, local_coords);
                if (interp.fn(forest, cell_idx, coef, interp.userdata) != MR_SUCCESS) {
                    return MR_FAILURE;
                }
            }
        }
    }

    return MR_SUCCESS;
}

#if 0
static mr_int descend_to_finest_neighbor(mr_ocforest *forest, mr_int node_idx, mr_direction dir[MR_ADJACENCY_VERTEX]) {
    if (node_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    while (!(node->flags & MR_OCTREE_NODE_FLAG_LEAF)) {
        node_idx = node->first_child + mr_direction_to_local_idx(dir);
        node = mr_ocforest_get_node(forest, node_idx);
    }

    return node_idx;
}

static int find_ghost_cell_stencil(mr_ocforest *forest, mr_int node_idx, mr_direction dir[MR_ADJACENCY_VERTEX], mr_int stencil[Q_STENCIL_DIM]) {
    mr_int full_stencil[Q_STENCIL_DIM + 1] = { MR_INVALID_INDEX, node_idx, MR_INVALID_INDEX, MR_INVALID_INDEX };
    mr_int stencil_shift = 0;

    mr_direction reflected_dir[] = {
        mr_direction_reflect(dir[0]),
        mr_direction_reflect(dir[1]),
        mr_direction_reflect(dir[2]),
    };

    full_stencil[0] = mr_octree_find_vertex_neighbor(forest, node_idx, dir);
    full_stencil[2] = mr_octree_find_vertex_neighbor(forest, node_idx, reflected_dir);

    full_stencil[0] = descend_to_finest_neighbor(forest, full_stencil[0], reflected_dir);
    full_stencil[2] = descend_to_finest_neighbor(forest, full_stencil[2], dir);

    if (full_stencil[0] == MR_INVALID_INDEX && full_stencil[2] == MR_INVALID_INDEX) {
        return MR_FAILURE;
    }

    if (full_stencil[0] == MR_INVALID_INDEX) {
        full_stencil[3] = mr_octree_find_vertex_neighbor(forest, full_stencil[2], reflected_dir);
        if (full_stencil[3] == MR_INVALID_INDEX) {
            return MR_FAILURE;
        }

        full_stencil[3] = descend_to_finest_neighbor(forest, full_stencil[3], dir);
        stencil_shift = 1;
    } else if (full_stencil[2] == MR_INVALID_INDEX) {
        full_stencil[3] = mr_octree_find_vertex_neighbor(forest, full_stencil[0], dir);
        if (full_stencil[3] == MR_INVALID_INDEX) {
            return MR_FAILURE;
        }

        full_stencil[3] = descend_to_finest_neighbor(forest, full_stencil[3], reflected_dir);
        stencil_shift = 1;
    }
    
    stencil[0] = full_stencil[MR_MOD(stencil_shift, Q_STENCIL_DIM)];
    stencil[1] = full_stencil[MR_MOD(stencil_shift + 1, Q_STENCIL_DIM)];
    stencil[2] = full_stencil[MR_MOD(stencil_shift + 2, Q_STENCIL_DIM)];

    return MR_SUCCESS;
}
#endif

int mr_fvm_calculate_ghost_cell(mr_ocforest *forest, mr_int node_idx, mr_direction dir[MR_ADJACENCY_VERTEX], mr_fvm_interpolation_cb interp) {
#if 0
    if (!forest || node_idx == MR_INVALID_INDEX || !dir || !interp.fn) {
        return MR_FAILURE;
    }

    mr_int stencil[Q_STENCIL_DIM] = { 0 };
    if (find_ghost_cell_stencil(forest, node_idx, dir, stencil) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    mr_float ghost_cell_center[] = { node->x, node->y, node->z };
    for (size_t i = 0; i < MR_ADJACENCY_VERTEX; ++i) {
        mr_axis axis = mr_direction_get_axis(dir[i]);
        mr_float sign_mul = mr_direction_get_sign_mul(dir[i]);

        ghost_cell_center[axis] += sign_mul * node->dim / 4.0f;
    }

    mr_octree_node *first_node = mr_ocforest_get_node(forest, stencil[0]);
    mr_float dist = mr_norm2(
        ghost_cell_center[MR_AXIS_X] - first_node->x,
        ghost_cell_center[MR_AXIS_Y] - first_node->y,
        ghost_cell_center[MR_AXIS_Z] - first_node->z
    );

    mr_float data_points[Q_STENCIL_DIM] = { 0.0f };
    for (size_t i = 1; i < Q_STENCIL_DIM; ++i) {
        node = mr_ocforest_get_node(forest, stencil[i]);
        data_points[i] = mr_norm2(node->x - first_node->x, node->y - first_node->y, node->z - first_node->z);
    }

    for (size_t i = 0; i < Q_STENCIL_DIM; ++i) {
        mr_float coef = lagrange_coef(dist, data_points, i);
        interp.fn(forest, stencil[i], coef, interp.userdata);
    }

    return MR_SUCCESS;
#else
    return MR_FAILURE;
#endif
}