#include <stdio.h>
#include <string.h>

#include "maniray/compute/fvm/interpolation.h"
#include "maniray/utils/misc.h"

// Q for triQuadratic interpolation
#define Q_STENCIL_DIM 3
#define Q_STENCIL_DIM_SQR (Q_STENCIL_DIM * Q_STENCIL_DIM)

#define WRAP_CNT(cnt, min, max) ((cnt) <= (max) ? (cnt) : (max) - (cnt))

// For loop which goes from 0 to max, then from -1 to min. 0 must be included in the [min; max] interval.
#define BIDIR_FOR(type, cnt, min, max) \
    for (type _pos_ ## cnt = 0, cnt = 0; cnt >= (min); ++_pos_ ## cnt, cnt = WRAP_CNT(_pos_ ## cnt, min, max))

static mr_uint flatten_idx(mr_int first_local_idx[MR_NB_AXES], mr_int local_idx[MR_NB_AXES]) {
    mr_uint pos_i = local_idx[MR_AXIS_X] - first_local_idx[MR_AXIS_X];
    mr_uint pos_j = local_idx[MR_AXIS_Y] - first_local_idx[MR_AXIS_Y];
    mr_uint pos_k = local_idx[MR_AXIS_Z] - first_local_idx[MR_AXIS_Z];

    return pos_i * Q_STENCIL_DIM_SQR + pos_j * Q_STENCIL_DIM + pos_k;
}

static bool is_node_visited(mr_bitfield visited_nodes, mr_int first_local_idx[MR_NB_AXES], mr_int local_idx[MR_NB_AXES]) {
    mr_uint flat_idx = flatten_idx(first_local_idx, local_idx);

    return (bool)(visited_nodes >> flat_idx & 1);
}

static mr_int get_next_node_index(mr_ocforest *forest, mr_int first_idx, mr_int prev_idx, mr_int offset, mr_axis axis) {
    if (offset == 0) {
        return first_idx;
    }

    mr_sign sign = offset < 0 ? MR_SIGN_MINUS : MR_SIGN_PLUS;
    mr_direction dir = mr_direction_create(axis, sign);

    if (offset == -1 || offset == 1) {
        return mr_octree_find_face_neighbor(forest, first_idx, dir);
    }

    return mr_octree_find_face_neighbor(forest, prev_idx, dir);
}

static mr_float lagrange_coef(mr_float x, mr_int first_idx, mr_int idx) {
    mr_float res = 1.0f;
    for (mr_int i = 0; i < Q_STENCIL_DIM; ++i) {
        res *= (x - (mr_float)(first_idx + i)) / (mr_float)(idx - first_idx - i);
    }

    return res;
}

static mr_float calculate_coef(
    const mr_float p[MR_NB_AXES],
    mr_octree_node *node,
    mr_int first_local_idx[MR_NB_AXES],
    mr_int local_idx[MR_NB_AXES]
) {
    if (!p || !node || !first_local_idx) {
        return 0.0f;
    }

    mr_float local_coords[] = {
        (p[MR_AXIS_X] - node->x) / node->dim,
        (p[MR_AXIS_Y] - node->y) / node->dim,
        (p[MR_AXIS_Z] - node->z) / node->dim,
    };

    mr_float coef = 1.0f;
    for (mr_int i = 0; i < Q_STENCIL_DIM; ++i) {
        coef *= lagrange_coef(local_coords[i], first_local_idx[i], local_idx[i]);
    }

    return coef;
}

static int q_stencil_apply(
    mr_ocforest *forest,
    const mr_float p[MR_NB_AXES],
    mr_int node_idx,
    mr_int first_local_idx[MR_NB_AXES],
    mr_fvm_interpolation_cb cb
) {
    if (node_idx == MR_INVALID_INDEX) {
        return MR_FAILURE;
    }

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    mr_bitfield visited_nodes = 0;

    mr_int x_node_idx = MR_INVALID_INDEX;
    mr_int y_node_idx = MR_INVALID_INDEX;
    mr_int z_node_idx = MR_INVALID_INDEX;
    BIDIR_FOR(mr_int, i, first_local_idx[0], first_local_idx[0] + 2) {
        x_node_idx = get_next_node_index(forest, node_idx, x_node_idx, i, MR_AXIS_X);
        if (x_node_idx == MR_INVALID_INDEX) {
            return MR_FAILURE;
        }

        BIDIR_FOR(mr_int, j, first_local_idx[1], first_local_idx[1] + 2) {
            y_node_idx = get_next_node_index(forest, x_node_idx, y_node_idx, j, MR_AXIS_Y);
            if (y_node_idx == MR_INVALID_INDEX) {
                return MR_FAILURE;
            }

            BIDIR_FOR(mr_int, k, first_local_idx[2], first_local_idx[2] + 2) {
                mr_int local_idx[] = { i, j, k };
                if (is_node_visited(visited_nodes, first_local_idx, local_idx)) {
                    continue;
                }

                z_node_idx = get_next_node_index(forest, y_node_idx, z_node_idx, k, MR_AXIS_Z);
                if (z_node_idx == MR_INVALID_INDEX) {
                    return MR_FAILURE;
                }

                mr_float coef = calculate_coef(p, node, first_local_idx, local_idx);
                cb.fn(forest, z_node_idx, local_idx, coef, cb.userdata);

                visited_nodes |= 1 << flatten_idx(first_local_idx, local_idx);
            }
        }
    }

    return MR_SUCCESS;
}

typedef struct stencil_available_userdata {
    mr_octree_node *host_node;
    bool is_available;
} stencil_available_userdata;

static int check_q_stencil_available(mr_ocforest *forest, mr_int node_idx, mr_int local_idx[MR_NB_AXES], mr_float coef, void *userdata) {
    MR_UNUSED(local_idx);
    MR_UNUSED(coef);

    stencil_available_userdata *data = userdata;
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);

    if (!(node->flags & MR_OCTREE_NODE_FLAG_ACTIVE) || !(node->flags & MR_OCTREE_NODE_FLAG_LEAF) || data->host_node->level != node->level) {
        data->is_available = false;
    }

    return MR_SUCCESS;
}

static int find_available_q_stencil(
    mr_ocforest *forest,
    mr_int node_idx,
    mr_octree_node *node,
    mr_int local_idx[MR_NB_AXES],
    mr_int res_local_idx[MR_NB_AXES]
) {
    stencil_available_userdata data = { node, true };
    int res = q_stencil_apply(forest, NULL, node_idx, local_idx, mr_fvm_interpolation_cb_create(check_q_stencil_available, &data));

    if (data.is_available && res == MR_SUCCESS) {
        memcpy(res_local_idx, local_idx, MR_NB_AXES * sizeof(mr_int));

        return MR_SUCCESS;
    }

    return MR_FAILURE;
}

static int find_available_center_q_stencil(mr_ocforest *forest, mr_int node_idx, mr_octree_node *node, mr_int res_local_idx[MR_NB_AXES]) {
    mr_int local_idx[MR_NB_AXES] = { -1, -1, -1 };

    return find_available_q_stencil(forest, node_idx, node, local_idx, res_local_idx);
}

static int find_available_face_q_stencil(mr_ocforest *forest, mr_int node_idx, mr_octree_node *node, mr_int res_local_idx[MR_NB_AXES]) {
    for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
        mr_int local_idx[] = { -1, -1, -1 };
        local_idx[mr_direction_get_axis(dir)] += mr_direction_get_sign_mul(dir);

        if (find_available_q_stencil(forest, node_idx, node, local_idx, res_local_idx) == MR_SUCCESS) {
            return MR_SUCCESS;
        }
    }

    return MR_FAILURE;
}

static int find_available_edge_q_stencil(mr_ocforest *forest, mr_int node_idx, mr_octree_node *node, mr_int res_local_idx[MR_NB_AXES]) {
    for (mr_direction dir_0 = MR_DIRECTION_MI_X; dir_0 <= MR_DIRECTION_PL_Y; ++dir_0) {
        for (mr_direction dir_1 = dir_0 - mr_direction_get_sign(dir_0) + 2; dir_1 <= MR_DIRECTION_PL_Z; ++dir_1) {
            mr_int local_idx[] = { -1, -1, -1 };
            local_idx[mr_direction_get_axis(dir_0)] += mr_direction_get_sign_mul(dir_0);
            local_idx[mr_direction_get_axis(dir_1)] += mr_direction_get_sign_mul(dir_1);

            if (find_available_q_stencil(forest, node_idx, node, local_idx, res_local_idx) == MR_SUCCESS) {
                return MR_SUCCESS;
            }
        }
    }

    return MR_FAILURE;
}

static int find_available_vertex_q_stencil(mr_ocforest *forest, mr_int node_idx, mr_octree_node *node, mr_int res_local_idx[MR_NB_AXES]) {
#define NB_VERTICES 8
    for (mr_uint i = 0; i < NB_VERTICES; ++i) {
        mr_sign signs[] = {
            i & 1,
            i >> 1 & 1,
            i >> 2 & 1,
        };

        mr_int local_idx[] = { -1, -1, -1 };
        for (mr_int i = 0; i < MR_NB_AXES; ++i) {
            local_idx[i] += mr_sign_to_mul(signs[i]);
        }

        if (find_available_q_stencil(forest, node_idx, node, local_idx, res_local_idx) == MR_SUCCESS) {
            return MR_SUCCESS;
        }
    }

    return MR_FAILURE;
}

static int find_q_stencil(
    mr_ocforest *forest,
    mr_int node_idx,
    mr_int res_local_idx[MR_NB_AXES]
) {
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    if (!node) {
        return MR_FAILURE;
    }

    if (find_available_center_q_stencil(forest, node_idx, node, res_local_idx) == MR_SUCCESS) {
        return MR_SUCCESS;
    }

    if (find_available_face_q_stencil(forest, node_idx, node, res_local_idx) == MR_SUCCESS) {
        return MR_SUCCESS;
    }

    if (find_available_edge_q_stencil(forest, node_idx, node, res_local_idx) == MR_SUCCESS) {
        return MR_SUCCESS;
    }

    if (find_available_vertex_q_stencil(forest, node_idx, node, res_local_idx) == MR_SUCCESS) {
        return MR_SUCCESS;
    }

    return MR_FAILURE;
}

bool mr_fvm_point_has_interpolation_stencil(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES]) {
    if (!forest) {
        return false;
    }

    mr_int host_node_idx = mr_octree_locate_point(forest, octree_idx, p);
    mr_int local_idx[MR_NB_AXES] = { 0 };

    return find_q_stencil(forest, host_node_idx, local_idx) == MR_SUCCESS;
}

int mr_fvm_perform_interpolation(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES], mr_fvm_interpolation_cb interp) {
    if (!forest || !p || !interp.fn) {
        return MR_FAILURE;
    }

    mr_int local_idx[MR_NB_AXES] = { 0 };
    mr_int host_node_idx = mr_octree_locate_point(forest, octree_idx, p);
    if (find_q_stencil(forest, host_node_idx, local_idx) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    return q_stencil_apply(forest, p, host_node_idx, local_idx, interp);
}