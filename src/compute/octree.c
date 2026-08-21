#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <tgmath.h>

#include "maniray/compute/math.h"
#include "maniray/compute/octree.h"
#include "maniray/utils/xmalloc.h"
#include "maniray/utils/misc.h"

static void calculate_cell_center(mr_int local_idx, mr_float node_dim, mr_float node_center[], mr_float cell_center[]) {
    mr_int i = local_idx & 0x3;
    mr_int j = local_idx >> 2 & 0x3;
    mr_int k = local_idx >> 4 & 0x3;

    mr_float cell_dim = node_dim / MR_OCTREE_NODE_BLOCK_DIM;

    cell_center[0] = node_center[0] + (i - 1.5) * cell_dim;
    cell_center[1] = node_center[1] + (j - 1.5) * cell_dim;
    cell_center[2] = node_center[2] + (k - 1.5) * cell_dim;
}

mr_ocforest *mr_ocforest_create(
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots,
    size_t cell_extra_fields[],
    size_t nb_cell_extra_fields
) {
    mr_ocforest *forest = xmalloc(sizeof(mr_ocforest));

    forest->manifold = manifold;
    forest->nb_roots = nb_roots;
    forest->roots = xmalloc(nb_roots * sizeof(mr_octree_root));

    forest->nodes = mr_mem_pool_create(
        (size_t[]) { sizeof(mr_octree_node), sizeof(mr_octree_node_connection) },
        MR_OCTREE_NODE_NB_MAIN_FIELDS
    );

    size_t nb_fields = MR_OCTREE_CELL_NB_MAIN_FIELDS + nb_cell_extra_fields;
    size_t *field_sizes = xmalloc(nb_fields * sizeof(size_t));

    field_sizes[0] = sizeof(mr_octree_cell);
    memcpy(field_sizes + MR_OCTREE_CELL_NB_MAIN_FIELDS, cell_extra_fields, nb_cell_extra_fields * sizeof(size_t));

    forest->cells = mr_mem_pool_create(field_sizes, nb_fields);
    free(field_sizes);

    mr_index index = mr_mem_pool_alloc_many(forest->nodes, nb_roots);
    assert(index <= INT32_MAX && index + MR_OCTREE_NB_CHILDREN <= INT32_MAX);

    mr_int node_idx = (mr_int)index;
    if (node_idx == MR_INVALID_INDEX) {
        goto failure;
    }

    for (mr_int i = 0; (size_t)i < nb_roots; ++i) {
        forest->roots[i] = (mr_octree_root) {
            .flags = roots[i].flags,
            .node_idx = node_idx + i,
        };

        mr_octree_node *node = mr_ocforest_get_node(forest, node_idx + i);
        *node = (mr_octree_node) {
            .flags = MR_OCTREE_NODE_FLAG_ACTIVE | MR_OCTREE_NODE_FLAG_LEAF,
            .level = 0,
            .chart_idx = roots[i].chart_idx,
            .x = roots[i].x,
            .y = roots[i].y,
            .z = roots[i].z,
            .dim = roots[i].dim,
            .root = i,
            .parent = MR_INVALID_INDEX,
            .first_child = MR_INVALID_INDEX,
            .value = 0.0f,
            .gradient = { 0.0f },
        };

        mr_octree_node_connection *conn = mr_ocforest_get_node_connection(forest, node_idx + i);
        *conn = (mr_octree_node_connection) {
            .local_idx = MR_INVALID_INDEX,
            .external_neighbors = { MR_INVALID_INDEX, MR_INVALID_INDEX, MR_INVALID_INDEX },
        };

        index = mr_mem_pool_alloc_many(forest->cells, MR_OCTREE_NB_CELLS_IN_BLOCK);
        assert(index <= INT32_MAX && index + MR_OCTREE_NB_CELLS_IN_BLOCK <= INT32_MAX);

        mr_int first_cell_idx = (mr_int)index;
        if (first_cell_idx == MR_INVALID_INDEX) {
            goto failure;
        }
        node->first_child = first_cell_idx;

        for (mr_int local_cell_idx = 0; local_cell_idx < MR_OCTREE_NB_CELLS_IN_BLOCK; ++local_cell_idx) {
            mr_octree_cell *cell = mr_ocforest_get_cell(forest, first_cell_idx + local_cell_idx);

            mr_float center[MR_NB_AXES] = { 0.0f };
            calculate_cell_center(local_cell_idx, node->dim, (mr_float[]) { node->x, node->y, node->z }, center);

            *cell = (mr_octree_cell) {
                .parent = node_idx + i,
                .chart_idx = node->chart_idx,
                .x = center[0],
                .y = center[1],
                .z = center[2],
                .dim = node->dim / MR_OCTREE_NODE_BLOCK_DIM,
            };
        }
    }

    return forest;

failure:
    mr_mem_pool_destroy(forest->cells);
    mr_mem_pool_destroy(forest->nodes);
    free(forest->roots);
    free(forest);

    return NULL;
}

void mr_ocforest_destroy(mr_ocforest *forest) {
    if (!forest) {
        return;
    }

    mr_mem_pool_destroy(forest->cells);
    mr_mem_pool_destroy(forest->nodes);
    free(forest->roots);
    free(forest);
}

size_t mr_ocforest_nb_nodes_upper_bound(mr_ocforest *forest) {
    return forest ? mr_mem_pool_len_bound(forest->nodes) : 0;
}

size_t mr_ocforest_nb_cells_upper_bound(mr_ocforest *forest) {
    return forest ? mr_mem_pool_len_bound(forest->cells) : 0;
}

static int leaves_counter(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    MR_UNUSED(forest);
    MR_UNUSED(node_idx);

    ++(*(size_t *)userdata);

    return MR_SUCCESS;
}

size_t mr_ocforest_count_leaves(mr_ocforest *forest) {
    assert(forest);

    size_t nb_leaves = 0;
    for (mr_index octree_idx = 0; (size_t)octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_leaves_apply(forest, octree_idx, mr_octree_apply_cb_create(leaves_counter, &nb_leaves));
    }

    return nb_leaves;
}

size_t mr_ocforest_count_cells(mr_ocforest *forest) {
    return mr_ocforest_count_leaves(forest) * MR_OCTREE_NB_CELLS_IN_BLOCK;
}

mr_octree_node *mr_ocforest_get_node(mr_ocforest *forest, mr_int idx) {
    assert(forest);

    return mr_mem_pool_ptr(forest->nodes, MR_OCTREE_NODE_FIELD, idx);
}

mr_octree_node *mr_ocforest_get_node_array(mr_ocforest *forest) {
    assert(forest);

    return mr_mem_pool_array_ptr(forest->nodes, MR_OCTREE_NODE_FIELD);
}

mr_octree_node_connection *mr_ocforest_get_node_connection(mr_ocforest *forest, mr_int idx) {
    assert(forest);

    return mr_mem_pool_ptr(forest->nodes, MR_OCTREE_NODE_CONNECTION_FIELD, idx);
}

mr_octree_node_connection *mr_ocforest_get_node_connection_array(mr_ocforest *forest) {
    assert(forest);

    return mr_mem_pool_array_ptr(forest->nodes, MR_OCTREE_NODE_CONNECTION_FIELD);
}

mr_octree_cell *mr_ocforest_get_cell(mr_ocforest *forest, mr_int idx) {
    assert(forest);

    return mr_mem_pool_ptr(forest->cells, MR_OCTREE_CELL_FIELD, idx);
}

mr_octree_cell *mr_ocforest_get_cell_array(mr_ocforest *forest) {
    assert(forest);

    return mr_mem_pool_array_ptr(forest->cells, MR_OCTREE_CELL_FIELD);
}

void *mr_ocforest_get_cell_extra(mr_ocforest *forest, mr_int idx, mr_int field) {
    assert(forest);

    return mr_mem_pool_ptr(forest->cells, MR_OCTREE_CELL_NB_MAIN_FIELDS + field, idx);
}

void *mr_ocforest_get_cell_extra_array(mr_ocforest *forest, mr_int field) {
    assert(forest);

    return mr_mem_pool_array_ptr(forest->cells, MR_OCTREE_CELL_NB_MAIN_FIELDS + field);
}

static void get_root_corner(mr_octree_node *root, mr_float coords[MR_NB_AXES]) {
    mr_float hdim = root->dim / 2.0f;

    coords[MR_AXIS_X] = root->x - hdim;
    coords[MR_AXIS_Y] = root->y - hdim;
    coords[MR_AXIS_Z] = root->z - hdim;
}

static void get_node_int_coords(mr_ocforest *forest, mr_octree_node *node, mr_octree_cell *cell, mr_int coords[MR_NB_AXES]) {
    mr_octree_node *root = mr_ocforest_get_node(forest, forest->roots[node->root].node_idx);

    mr_float root_coords[MR_NB_AXES] = { 0.0f };
    get_root_corner(root, root_coords);    

    mr_float hdim = cell->dim / 2.0f;
    mr_float scale = (mr_float)(1 << (MR_OCTREE_MAX_LEVEL + 2)) / root->dim;

    coords[MR_AXIS_X] = (mr_int)llround((cell->x - root_coords[MR_AXIS_X] - hdim) * scale);
    coords[MR_AXIS_Y] = (mr_int)llround((cell->y - root_coords[MR_AXIS_Y] - hdim) * scale);
    coords[MR_AXIS_Z] = (mr_int)llround((cell->z - root_coords[MR_AXIS_Z] - hdim) * scale);

    mr_int max_value = (1 << (MR_OCTREE_MAX_LEVEL + 2)) - 1;
    coords[MR_AXIS_X] = MR_CLAMP(coords[MR_AXIS_X], 0, max_value);
    coords[MR_AXIS_Y] = MR_CLAMP(coords[MR_AXIS_Y], 0, max_value);
    coords[MR_AXIS_Z] = MR_CLAMP(coords[MR_AXIS_Z], 0, max_value);
}

static mr_int interleave_bits(mr_int ix, mr_int iy, mr_int iz, size_t bits) {
    mr_int code = 0;
    for (size_t b = 0; b < bits; b++) {
        code |= ((mr_int)((ix >> b) & 1) << (3 * b));
        code |= ((mr_int)((iy >> b) & 1) << (3 * b + 1));
        code |= ((mr_int)((iz >> b) & 1) << (3 * b + 2));
    }

    return code;
}

mr_int mr_ocforest_get_code(mr_ocforest *forest, mr_int cell_idx) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);
    assert(1 << (MR_INT_NB_BITS - MR_NB_AXES * (MR_OCTREE_MAX_LEVEL + 2) - 1) > forest->nb_roots);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);

    mr_int coords[MR_NB_AXES] = { 0 };
    get_node_int_coords(forest, node, cell, coords);

    mr_int morton_code = interleave_bits(coords[MR_AXIS_X], coords[MR_AXIS_Y], coords[MR_AXIS_Z], MR_OCTREE_MAX_LEVEL);
    mr_int root_code = node->root << MR_NB_AXES * MR_OCTREE_MAX_LEVEL;

    return morton_code | root_code;
}

static mr_int extract_local_idx(mr_int code, mr_int level) {
    return code >> MR_NB_AXES * (MR_OCTREE_MAX_LEVEL + 2 - level) & 0x7;
}

mr_int mr_ocforest_find_cell_with_code(mr_ocforest *forest, mr_int code) {
    assert(forest);
    assert(code != MR_INVALID_INDEX);

    mr_int root_idx = code >> MR_NB_AXES * MR_OCTREE_MAX_LEVEL;

    mr_int node_idx = forest->roots[root_idx].node_idx;
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);

    for (mr_int level = 1; level <= MR_OCTREE_MAX_LEVEL; ++level) {
        if (node->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            return node->first_child
                + (extract_local_idx(code, level + 1) << MR_NB_AXES)
                + extract_local_idx(code, level + 2);
        }

        node_idx = node->first_child + extract_local_idx(code, level);
        node = mr_ocforest_get_node(forest, node_idx);
    }

    return MR_INVALID_INDEX;
}

static int mr_octree_leaves_apply_ext(mr_ocforest *forest, mr_index octree_idx, mr_octree_apply_cb apply, bool recursive) {
    assert(forest);
    assert((size_t)octree_idx < forest->nb_roots);

    mr_int root_idx = forest->roots[octree_idx].node_idx;
    mr_octree_node *root = mr_ocforest_get_node(forest, root_idx);
    if (root->flags & MR_OCTREE_NODE_FLAG_LEAF) {
        apply.fn(forest, root_idx, apply.userdata);

        if (!recursive) {
            return MR_SUCCESS;
        }
    }

    mr_int parent_idx = root_idx;
    mr_octree_node *parent = NULL;

    mr_int child_idx = 0;
    mr_octree_node *child = NULL;
    mr_uint saved_flags = 0;

    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };
    mr_int level = 0;
    while (true) {
        parent = mr_ocforest_get_node(forest, parent_idx);
        if (path[level] == MR_OCTREE_NB_CHILDREN) {
            if (--level < 0) {
                break;
            }

            ++path[level];
            parent_idx = parent->parent;

            continue;
        }

        child_idx = parent->first_child + path[level];
        saved_flags = mr_ocforest_get_node(forest, child_idx)->flags;
        if (saved_flags & MR_OCTREE_NODE_FLAG_LEAF) {
            if (apply.fn(forest, child_idx, apply.userdata) == MR_FAILURE) {
                return MR_FAILURE;
            }
        }

        child = mr_ocforest_get_node(forest, child_idx);
        if ((saved_flags & MR_OCTREE_NODE_FLAG_LEAF && !recursive) || child->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            ++path[level];
        } else {
            parent_idx = child_idx;
            path[++level] = 0;
        }
    }

    return MR_SUCCESS;
}

int mr_octree_leaves_apply(mr_ocforest *forest, mr_index octree_idx, mr_octree_apply_cb apply) {
    return mr_octree_leaves_apply_ext(forest, octree_idx, apply, false);
}

static int apply_cells_in_leaf(mr_ocforest *forest, mr_int idx, void *userdata) {
    mr_octree_apply_cb *cb = userdata;
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    for (mr_int local_cell_idx = 0; local_cell_idx < MR_OCTREE_NB_CELLS_IN_BLOCK; ++local_cell_idx) {
        if (cb->fn(forest, node->first_child + local_cell_idx, cb->userdata) != MR_SUCCESS) {
            return MR_FAILURE;
        }
    }

    return MR_SUCCESS;
}

int mr_octree_cells_apply(mr_ocforest *forest, mr_index octree_idx, mr_octree_apply_cb apply) {
    return mr_octree_leaves_apply(forest, octree_idx, mr_octree_apply_cb_create(apply_cells_in_leaf, &apply));
}

void mr_octree_periodic_wrap(mr_ocforest *forest, mr_index octree_idx, mr_float p[MR_NB_AXES]) {
    assert(forest);
    assert(p);
    assert((size_t)octree_idx < forest->nb_roots);

    mr_octree_root *root = &forest->roots[octree_idx];
    mr_octree_node *root_node = mr_ocforest_get_node(forest, root->node_idx);

    mr_float hdim = root_node->dim / 2.0f;

    if (root->flags & MR_OCTREE_FLAG_PERIODIC_X) {
        p[0] = mr_wrap(p[0], root_node->x - hdim, root_node->x + hdim);
    }

    if (root->flags & MR_OCTREE_FLAG_PERIODIC_Y) {
        p[1] = mr_wrap(p[1], root_node->y - hdim, root_node->y + hdim);
    }

    if (root->flags & MR_OCTREE_FLAG_PERIODIC_Z) {
        p[2] = mr_wrap(p[2], root_node->z - hdim, root_node->z + hdim);
    }
}

mr_int mr_octree_locate_point_in_leaf(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES]) {
    assert(forest);
    assert(p);
    assert((size_t)octree_idx < forest->nb_roots);

    mr_int idx = forest->roots[octree_idx].node_idx;
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    mr_float wp[3] = { p[0], p[1], p[2] };
    mr_octree_periodic_wrap(forest, octree_idx, wp);

    mr_int child_idx = 0;
    for (size_t level = 0; level <= MR_OCTREE_MAX_LEVEL; ++level) {
        node = mr_ocforest_get_node(forest, idx);
        if (node->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            break;
        }

        child_idx = (mr_int)(wp[0] > node->x)
            | (mr_int)(wp[1] > node->y) << 1
            | (mr_int)(wp[2] > node->z) << 2;
        idx = node->first_child + child_idx;
    }

    if (mr_norm_inf(wp[0] - node->x, wp[1] - node->y, wp[2] - node->z) > node->dim / 2.0f) {
        return MR_INVALID_INDEX;
    }

    return idx;
}

mr_int mr_octree_locate_point_in_cell(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES]) {
    assert(forest);
    assert(p);
    assert((size_t)octree_idx < forest->nb_roots);

    mr_float wp[3] = { p[0], p[1], p[2] };
    mr_octree_periodic_wrap(forest, octree_idx, wp);

    mr_int node_idx = mr_octree_locate_point_in_leaf(forest, octree_idx, wp);
    if (node_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);

    mr_float hdim = node->dim / 2.0f;
    mr_float cell_dim = node->dim / MR_OCTREE_NODE_BLOCK_DIM;

    mr_int local_coords[] = {
        (mr_int)floor((wp[0] - node->x + hdim) / cell_dim),
        (mr_int)floor((wp[1] - node->y + hdim) / cell_dim),
        (mr_int)floor((wp[2] - node->z + hdim) / cell_dim),
    };

    local_coords[0] = MR_CLAMP(local_coords[0], 0, 3);
    local_coords[1] = MR_CLAMP(local_coords[1], 0, 3);
    local_coords[2] = MR_CLAMP(local_coords[2], 0, 3);

    mr_int local_idx = (local_coords[2] << 4) + (local_coords[1] << 2) + local_coords[0];

    return node->first_child + local_idx;
}

static int insert_children(mr_ocforest *forest, mr_int parent_idx, void *userdata);

static mr_int get_cached_neighbor(mr_ocforest *forest, mr_int node_idx, mr_direction dir, bool *might_change) {
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    mr_octree_node_connection *conn = mr_ocforest_get_node_connection(forest, node_idx);

    if (conn->local_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_axis axis = mr_direction_get_axis(dir);
    mr_sign sign = mr_direction_get_sign(dir);

    if ((conn->local_idx >> axis & 1) != sign) {
        mr_octree_node *parent = mr_ocforest_get_node(forest, node->parent);

        *might_change = false;
        return parent->first_child + (conn->local_idx ^ 1 << axis);
    }

    mr_int neighbor_idx = conn->external_neighbors[axis];
    if (neighbor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_octree_node *neighbor = mr_ocforest_get_node(forest, neighbor_idx);
    if (neighbor->level == node->level || neighbor->flags & MR_OCTREE_NODE_FLAG_LEAF) {
        *might_change = neighbor->level != node->level;
        return neighbor_idx;
    }

    return MR_INVALID_INDEX;
}

static void cache_neighbor(mr_ocforest *forest, mr_int node_idx, mr_int neighbor_idx, mr_direction dir) {
    mr_octree_node_connection *conn = mr_ocforest_get_node_connection(forest, node_idx);
    if (conn->local_idx == MR_INVALID_INDEX) {
        return;
    }

    mr_axis axis = mr_direction_get_axis(dir);
    mr_sign sign = mr_direction_get_sign(dir);

    if ((conn->local_idx >> axis & 1) != sign) {
        return;
    }

    conn->external_neighbors[axis] = neighbor_idx;
}

static mr_int get_local_idx(mr_ocforest *forest, mr_int node_idx) {
    mr_octree_node_connection *conn = mr_ocforest_get_node_connection(forest, node_idx);

    return conn->local_idx;
}

static bool is_node_adjacent(mr_int local_idx, mr_direction dir[], mr_adjacency adj) {
    for (size_t i = 0; i < adj; ++i) {
        mr_axis axis = mr_direction_get_axis(dir[i]);
        mr_sign sign = mr_direction_get_sign(dir[i]);

        if ((local_idx >> axis & 1) != sign) {
            return false;
        }
    }

    return true;
}

static bool is_root_adjacent(mr_ocforest *forest, mr_int root_idx, mr_direction dir[], mr_adjacency adj) {
    mr_octree_root *root = &forest->roots[root_idx];
    for (size_t i = 0; i < adj; ++i) {
        mr_axis axis = mr_direction_get_axis(dir[i]);

        if (!(root->flags & (1 << axis))) {
            return false;
        }
    }

    return true;
}

static mr_int local_idx_reflect_axes(mr_int local_idx, mr_direction dir[], mr_adjacency adj) {
    for (size_t i = 0; i < adj; ++i) {
        local_idx ^= 1 << mr_direction_get_axis(dir[i]);
    }

    return local_idx;
}

static mr_int find_ancestor(
    mr_ocforest *forest,
    mr_int node_idx,
    mr_direction dir[],
    mr_adjacency adj,
    mr_int path[MR_OCTREE_MAX_LEVEL],
    mr_uint *len
) {
    *len = 0;

    mr_int local_idx = MR_INVALID_INDEX;
    mr_int curr_idx = node_idx;
    mr_octree_node *node = NULL;
    do {
        node = mr_ocforest_get_node(forest, curr_idx);
        if (node->parent == MR_INVALID_INDEX) {
            if (is_root_adjacent(forest, node->root, dir, adj)) {
                break;
            }

            return MR_INVALID_INDEX;
        }

        local_idx = get_local_idx(forest, curr_idx);
        path[(*len)++] = local_idx;
        curr_idx = node->parent;
    } while (is_node_adjacent(local_idx, dir, adj));

    return curr_idx;
}

static void reflect_path(
    mr_direction dir[],
    mr_adjacency adj,
    mr_int path[MR_OCTREE_MAX_LEVEL],
    mr_uint len
) {
    for (mr_uint i = 0; i < len; ++i) {
        path[i] = local_idx_reflect_axes(path[i], dir, adj);
    }
}

static mr_int descend_path(
    mr_ocforest *forest,
    mr_int node_idx,
    mr_int path[MR_OCTREE_MAX_LEVEL],
    mr_uint len,
    bool do_refine,
    mr_int expected_level
) {
    mr_int curr_idx = node_idx;
    mr_octree_node *node = NULL;

    while (len > 0) {
        node = mr_ocforest_get_node(forest, curr_idx);
        if (do_refine && node->flags & MR_OCTREE_NODE_FLAG_LEAF && (mr_int)node->level < expected_level) {
            insert_children(forest, curr_idx, NULL);
            node = mr_ocforest_get_node(forest, curr_idx);
        }

        if (node->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            break;
        }

        curr_idx = node->first_child + path[--len];
    }

    return curr_idx;
}

static mr_int mr_octree_find_face_neighbor_node_with_refine(
    mr_ocforest *forest,
    mr_int node_idx,
    mr_direction dir,
    bool do_refine,
    mr_int expected_level
) {
    assert(forest);
    assert(node_idx != MR_INVALID_INDEX);

    bool cache_might_change = true;
    mr_int nidx = get_cached_neighbor(forest, node_idx, dir, &cache_might_change);
    if (nidx != MR_INVALID_INDEX && (!cache_might_change || !do_refine)) {
        return nidx;
    }

    mr_uint len = 0;
    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };

    mr_int ancestor_idx = find_ancestor(forest, node_idx, &dir, MR_ADJACENCY_FACE, path, &len);
    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    reflect_path(&dir, MR_ADJACENCY_FACE, path, len);
    nidx = descend_path(forest, ancestor_idx, path, len, do_refine, expected_level);
    cache_neighbor(forest, node_idx, nidx, dir);

    return nidx;
}

mr_int mr_octree_find_face_neighbor_node(mr_ocforest *forest, mr_int idx, mr_direction dir) {
    return mr_octree_find_face_neighbor_node_with_refine(forest, idx, dir, false, 0);
}

static bool is_balanced(mr_ocforest *forest, mr_int cell_idx, mr_int neighbor_idx, mr_direction dir) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_octree_node *neighbor_node = mr_ocforest_get_node(forest, neighbor_idx);

    return node->level - neighbor_node->level <= 1;
}

static mr_int read_cell_local_coord(mr_int local_idx, mr_axis axis) {
    return local_idx >> (axis * 2) & 0x3;
}

static mr_int write_cell_local_coord(mr_int local_idx, mr_axis axis, mr_int value) {
    local_idx &= ~(0x3 << axis * 2);
    local_idx |= (value & 0x3) << axis * 2;

    return local_idx;
}

static mr_int find_internal_neighbor_cell(mr_ocforest *forest, mr_int cell_idx, mr_direction dir) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);

    mr_int local_idx = cell_idx - node->first_child;
    mr_axis axis = mr_direction_get_axis(dir);
    mr_sign sign = mr_direction_get_sign(dir);

    mr_int local_coord = read_cell_local_coord(local_idx, axis);
    if (local_coord == (mr_int)sign * 0x3) {
        return MR_INVALID_INDEX;
    }

    local_coord += mr_sign_to_mul(sign);
    local_idx = write_cell_local_coord(local_idx, axis, local_coord);

    return node->first_child + local_idx;
}

static mr_int find_equal_size_neighbor_cell(mr_ocforest *forest, mr_int cell_idx, mr_int neighbor_idx, mr_direction dir) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_octree_node *neighbor_node = mr_ocforest_get_node(forest, neighbor_idx);
    
    if (node->level != neighbor_node->level || !(neighbor_node->flags & MR_OCTREE_NODE_FLAG_LEAF)) {
        return MR_INVALID_INDEX;
    }

    mr_axis axis = mr_direction_get_axis(dir);
    mr_int local_idx = (cell_idx - node->first_child) ^ (0x3 << axis * 2);

    return neighbor_node->first_child + local_idx;
}

static mr_int cell_idx_even_bits(mr_int cell_idx) {
    return (cell_idx >> 2 & 0x4) | (cell_idx >> 1 & 0x2) | (cell_idx & 0x1);
}

static mr_int cell_idx_odd_bits(mr_int cell_idx) {
    return (cell_idx >> 3 & 0x4) | (cell_idx >> 2 & 0x2) | (cell_idx >> 1 & 0x1);
}

// The local index of the cell is two interleaved local indices of octree nodes if the uniform grid of cells was given a tree structure
static mr_int interleave_local_idx(mr_int a, mr_int b) {
    return ((a & 0x4) << 3) | ((b & 0x4) << 2)
         | ((a & 0x2) << 2) | ((b & 0x2) << 1)
         | ((a & 0x1) << 1) | (b & 0x1);
}

static mr_int find_coarser_neighbor_cell(mr_ocforest *forest, mr_int cell_idx, mr_int neighbor_idx, mr_direction dir) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_octree_node *neighbor_node = mr_ocforest_get_node(forest, neighbor_idx);
    
    if (node->level - neighbor_node->level != 1 || !(neighbor_node->flags & MR_OCTREE_NODE_FLAG_LEAF)) {
        return MR_INVALID_INDEX;
    }

    mr_int node_local_idx = get_local_idx(forest, cell->parent);
    mr_int cell_outer_local_idx = cell_idx_odd_bits(cell_idx - node->first_child);

    mr_axis axis = mr_direction_get_axis(dir);
    mr_int local_idx = interleave_local_idx(node_local_idx ^ (1 << axis), cell_outer_local_idx ^ (1 << axis));

    return neighbor_node->first_child + local_idx;
}

static void find_finer_neighbor_cells(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_int neighbor_idx,
    mr_direction dir,
    mr_int neighbor_cells[4]
) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_octree_node *neighbor_node = mr_ocforest_get_node(forest, neighbor_idx);
    
    mr_axis axis = mr_direction_get_axis(dir);
    mr_sign sign = mr_direction_get_sign(dir);

    mr_int cell_local_idx = cell_idx - node->first_child;
    mr_int fine_neighbor_local_idx = cell_idx_odd_bits(cell_local_idx) ^ (1 << axis);

    mr_octree_node *fine_neighbor_node = mr_ocforest_get_node(forest, neighbor_node->first_child + fine_neighbor_local_idx);
    assert(fine_neighbor_node->first_child & MR_OCTREE_NODE_FLAG_LEAF);

    mr_int neighbor_cell_outer_idx = cell_idx_even_bits(cell_local_idx) ^ (1 << axis);
    size_t cnt = 0;
    for (mr_int neighbor_cell_inner_idx = 0; neighbor_cell_inner_idx < MR_OCTREE_NB_CHILDREN; ++neighbor_cell_inner_idx) {
        if ((neighbor_cell_inner_idx >> axis & 1) == sign) {
            continue;
        }

        neighbor_cells[cnt++] = fine_neighbor_node->first_child + interleave_local_idx(neighbor_cell_outer_idx, neighbor_cell_inner_idx);
    }
}

mr_octree_cell_neighbor mr_octree_find_face_neighbor_cells(mr_ocforest *forest, mr_int cell_idx, mr_direction dir) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

    mr_int internal_neighbor_idx = find_internal_neighbor_cell(forest, cell_idx, dir);
    if (internal_neighbor_idx != MR_INVALID_INDEX) {
        return (mr_octree_cell_neighbor) {
            .type = MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE,
            .node_idx = cell->parent,
            .neighbor_idx = internal_neighbor_idx,
        };
    }

    mr_int neighbor_idx = mr_octree_find_face_neighbor_node(forest, cell->parent, dir);

    if (neighbor_idx == MR_INVALID_INDEX) {
        return (mr_octree_cell_neighbor) {
            .type = MR_OCTREE_CELL_NEIGHBOR_NONE,
            .node_idx = MR_INVALID_INDEX,
        };
    }
    assert(is_balanced(forest, cell_idx, neighbor_idx, dir));

    mr_int equal_size_neighbor_idx = find_equal_size_neighbor_cell(forest, cell_idx, neighbor_idx, dir);
    if (equal_size_neighbor_idx != MR_INVALID_INDEX) {
        return (mr_octree_cell_neighbor) {
            .type = MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE,
            .node_idx = neighbor_idx,
            .neighbor_idx = equal_size_neighbor_idx,
        };
    }

    mr_int coarser_neighbor_idx = find_coarser_neighbor_cell(forest, cell_idx, neighbor_idx, dir);
    if (coarser_neighbor_idx != MR_INVALID_INDEX) {
        return (mr_octree_cell_neighbor) {
            .type = MR_OCTREE_CELL_NEIGHBOR_COARSER,
            .node_idx = neighbor_idx,
            .neighbor_idx = coarser_neighbor_idx,
        };
    }

    mr_octree_cell_neighbor ret = { .type = MR_OCTREE_CELL_NEIGHBOR_FINER, .node_idx = neighbor_idx, };
    find_finer_neighbor_cells(forest, cell_idx, neighbor_idx, dir, ret.neighbor_indices);
    
    return ret;
}

static bool cells_need_refine(mr_ocforest *forest, mr_int node_idx, mr_octree_cond_cb cond) {
    if (!cond.fn) {
        return true;
    }

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    for (mr_int i = 0; i < MR_OCTREE_NB_CELLS_IN_BLOCK; ++i) {
        mr_int cell_idx = node->first_child + i;
        if (cond.fn(forest, cell_idx, cond.userdata)) {
            return true;
        }
    }

    return false;
}

static int insert_children(mr_ocforest *forest, mr_int parent_idx, void *userdata) {
    mr_octree_cond_cb *cond_ud = userdata;
    if (cond_ud && !cells_need_refine(forest, parent_idx, *cond_ud)) {
        return MR_SUCCESS;
    }

    mr_octree_node *parent = mr_ocforest_get_node(forest, parent_idx);
    if (parent->level >= MR_OCTREE_MAX_LEVEL) {
        return MR_SUCCESS;
    }

    mr_index first_node_idx = mr_mem_pool_alloc_many(forest->nodes, MR_OCTREE_NB_CHILDREN);
    if (first_node_idx == MR_INVALID_INDEX) {
        return MR_FAILURE;
    }
    assert(first_node_idx <= INT32_MAX && first_node_idx + MR_OCTREE_NB_CHILDREN <= INT32_MAX);

    parent = mr_ocforest_get_node(forest, parent_idx);
    mr_int prev_cells = parent->first_child;

    parent->flags &= ~MR_OCTREE_NODE_FLAG_LEAF;
    parent->first_child = first_node_idx;

    mr_float qdim = parent->dim / 4.0f;
    for (size_t i = 0; i < MR_OCTREE_NB_CHILDREN; ++i) {
        mr_index first_cell_idx = mr_mem_pool_alloc_many(forest->cells, MR_OCTREE_NB_CELLS_IN_BLOCK);
        assert(first_cell_idx <= INT32_MAX && first_cell_idx + MR_OCTREE_NB_CELLS_IN_BLOCK <= INT32_MAX);
        if (first_cell_idx == MR_INVALID_INDEX) {
            return MR_FAILURE;
        }

        mr_octree_node *node = mr_ocforest_get_node(forest, first_node_idx + i);
        *node = (mr_octree_node) {
            .flags = MR_OCTREE_NODE_FLAG_UNBALANCED | MR_OCTREE_NODE_FLAG_LEAF | (parent->flags & MR_OCTREE_NODE_FLAG_ACTIVE),
            .level = parent->level + 1,
            .chart_idx = parent->chart_idx,
            .x = parent->x + (i & 1 ? qdim : -qdim),
            .y = parent->y + (i & 2 ? qdim : -qdim),
            .z = parent->z + (i & 4 ? qdim : -qdim),
            .dim = parent->dim / 2.0f,
            .root = parent->root,
            .parent = parent_idx,
            .first_child = first_cell_idx,
            .value = 0.0f,
            .gradient = { 0.0f },
        };

        mr_octree_node_connection *conn = mr_ocforest_get_node_connection(forest, first_node_idx + i);
        *conn = (mr_octree_node_connection) {
            .local_idx = i,
            .external_neighbors = { MR_INVALID_INDEX, MR_INVALID_INDEX, MR_INVALID_INDEX },
        };

        for (mr_int local_cell_idx = 0; local_cell_idx < MR_OCTREE_NB_CELLS_IN_BLOCK; ++local_cell_idx) {
            mr_octree_cell *cell = mr_ocforest_get_cell(forest, first_cell_idx + local_cell_idx);

            mr_float center[MR_NB_AXES] = { 0.0f };
            calculate_cell_center(local_cell_idx, node->dim, (mr_float[]) { node->x, node->y, node->z }, center);

            *cell = (mr_octree_cell) {
                .parent = first_node_idx + i,
                .chart_idx = node->chart_idx,
                .x = center[0],
                .y = center[1],
                .z = center[2],
                .dim = node->dim / MR_OCTREE_NODE_BLOCK_DIM,
            };
        }
    }

    mr_mem_pool_remove_many(forest->cells, prev_cells, MR_OCTREE_NB_CELLS_IN_BLOCK);

    return MR_SUCCESS;
}

void mr_octree_refine(mr_ocforest *forest, mr_index octree_idx, mr_octree_cond_cb cond, bool recursive) {
    mr_octree_leaves_apply_ext(forest, octree_idx, mr_octree_apply_cb_create(insert_children, &cond), recursive);
}

void mr_octree_refine_all(mr_ocforest *forest, mr_index octree_idx, size_t nb_repeats) {
    for (size_t i = 0; i < nb_repeats; ++i) {
        mr_octree_refine(forest, octree_idx, mr_octree_cond_cb_null(), false);
    }
}

static int balance_level(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    mr_uint level = *(mr_uint*)userdata;

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    if (node->flags & (MR_OCTREE_NODE_FLAG_UNBALANCED | MR_OCTREE_NODE_FLAG_LEAF) && node->level == level) {
        mr_int refine_level = (mr_int)node->level - 1;

        for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
            mr_octree_find_face_neighbor_node_with_refine(forest, node_idx, dir, true, refine_level);
        }

        node = mr_ocforest_get_node(forest, node_idx);
        node->flags &= ~MR_OCTREE_NODE_FLAG_UNBALANCED;
    }

    return MR_SUCCESS;
}

void mr_octree_balance(mr_ocforest *forest, mr_index octree_idx) {
    assert(forest);
    assert((size_t)octree_idx < forest->nb_roots);

    for (mr_uint i = MR_OCTREE_MAX_LEVEL; i > 0; --i) {
        mr_octree_leaves_apply(forest, octree_idx, mr_octree_apply_cb_create(balance_level, &i));
    }
}