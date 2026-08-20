#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <tgmath.h>

#include "maniray/compute/math.h"
#include "maniray/compute/octree.h"
#include "maniray/utils/xmalloc.h"
#include "maniray/utils/misc.h"

mr_ocforest *mr_ocforest_create(
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots,
    size_t extra_fields[],
    size_t nb_extra_fields
) {
    mr_ocforest *forest = xmalloc(sizeof(mr_ocforest));

    forest->manifold = manifold;
    forest->nb_roots = nb_roots;
    forest->roots = xmalloc(nb_roots * sizeof(mr_octree_root));

    size_t nb_fields = MR_OCTREE_NODE_NB_MAIN_FIELDS + nb_extra_fields;
    size_t *field_sizes = xmalloc(nb_fields * sizeof(size_t));

    field_sizes[0] = sizeof(mr_octree_node);
    field_sizes[1] = sizeof(mr_octree_node_connection);
    memcpy(field_sizes + MR_OCTREE_NODE_NB_MAIN_FIELDS, extra_fields, nb_extra_fields * sizeof(size_t));

    forest->nodes = mr_mem_pool_create(field_sizes, nb_fields);
    free(field_sizes);

    mr_int node_idx = (mr_int)mr_mem_pool_alloc_many(forest->nodes, nb_roots);
    if (node_idx == MR_INVALID_INDEX) {
        mr_mem_pool_destroy(forest->nodes);
        free(forest->roots);
        free(forest);

        return NULL;
    }

    for (size_t i = 0; i < nb_roots; ++i) {
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
            .root = (mr_int)i,
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
    }

    return forest;
}

void mr_ocforest_destroy(mr_ocforest *forest) {
    if (!forest) {
        return;
    }

    mr_mem_pool_destroy(forest->nodes);
    free(forest->roots);
    free(forest);
}

size_t mr_ocforest_size(mr_ocforest *forest) {
    return forest ? mr_mem_pool_len_bound(forest->nodes) : 0;
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
        mr_octree_leaves_apply(forest, octree_idx, mr_octree_cond_cb_null(), mr_octree_apply_cb_create(leaves_counter, &nb_leaves), false);
    }

    return nb_leaves;
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

void *mr_ocforest_get_extra(mr_ocforest *forest, mr_int idx, mr_int field) {
    assert(forest);

    return mr_mem_pool_ptr(forest->nodes, MR_OCTREE_NODE_NB_MAIN_FIELDS + field, idx);
}

void *mr_ocforest_get_extra_array(mr_ocforest *forest, mr_int field) {
    assert(forest);

    return mr_mem_pool_array_ptr(forest->nodes, MR_OCTREE_NODE_NB_MAIN_FIELDS + field);
}

static void get_root_corner(mr_octree_node *root, mr_float coords[MR_NB_AXES]) {
    mr_float hdim = root->dim / 2.0f;

    coords[MR_AXIS_X] = root->x - hdim;
    coords[MR_AXIS_Y] = root->y - hdim;
    coords[MR_AXIS_Z] = root->z - hdim;
}

static void get_node_int_coords(mr_ocforest *forest, mr_octree_node *node, mr_int coords[MR_NB_AXES]) {
    mr_octree_node *root = mr_ocforest_get_node(forest, forest->roots[node->root].node_idx);

    mr_float root_coords[MR_NB_AXES] = { 0.0f };
    get_root_corner(root, root_coords);    

    mr_float hdim = node->dim / 2.0f;
    mr_float scale = (mr_float)(1 << MR_OCTREE_MAX_LEVEL) / root->dim;

    coords[MR_AXIS_X] = (mr_int)llround((node->x - root_coords[MR_AXIS_X] - hdim) * scale);
    coords[MR_AXIS_Y] = (mr_int)llround((node->y - root_coords[MR_AXIS_Y] - hdim) * scale);
    coords[MR_AXIS_Z] = (mr_int)llround((node->z - root_coords[MR_AXIS_Z] - hdim) * scale);

    mr_int max_value = (1 << MR_OCTREE_MAX_LEVEL) - 1;
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

mr_int mr_ocforest_get_code(mr_ocforest *forest, mr_int idx) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);
    assert(1 << (MR_INT_NB_BITS - MR_NB_AXES * MR_OCTREE_MAX_LEVEL - 1) > forest->nb_roots);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    mr_int coords[MR_NB_AXES] = { 0 };
    get_node_int_coords(forest, node, coords);

    mr_int morton_code = interleave_bits(coords[MR_AXIS_X], coords[MR_AXIS_Y], coords[MR_AXIS_Z], MR_OCTREE_MAX_LEVEL);
    mr_int root_code = node->root << MR_NB_AXES * MR_OCTREE_MAX_LEVEL;

    return morton_code | root_code;
}

static mr_int extract_local_idx_from_code(mr_int code, mr_int level) {
    return code >> MR_NB_AXES * (MR_OCTREE_MAX_LEVEL - level) & 0x7;
}

mr_int mr_ocforest_find_node_with_code(mr_ocforest *forest, mr_int code) {
    assert(forest);
    assert(code != MR_INVALID_INDEX);

    mr_int root_idx = code >> MR_NB_AXES * MR_OCTREE_MAX_LEVEL;

    mr_int node_idx = forest->roots[root_idx].node_idx;
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);

    for (mr_int level = 1; level <= MR_OCTREE_MAX_LEVEL; ++level) {
        if (node->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            break;
        }

        node_idx = node->first_child + extract_local_idx_from_code(code, level);
        node = mr_ocforest_get_node(forest, node_idx);
    }

    return node_idx;
}

int mr_octree_leaves_apply(
    mr_ocforest *forest,
    mr_index octree_idx,
    mr_octree_cond_cb filter,
    mr_octree_apply_cb apply,
    bool recursive
) {
    assert(forest);
    assert((size_t)octree_idx < forest->nb_roots);

    mr_int root_idx = forest->roots[octree_idx].node_idx;
    if (filter.fn && !filter.fn(forest, root_idx, filter.userdata)) {
        return MR_SUCCESS;
    }

    mr_octree_node *root = mr_ocforest_get_node(forest, root_idx);
    if (root->first_child == MR_INVALID_INDEX) {
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
        if (filter.fn && !filter.fn(forest, child_idx, filter.userdata)) {
            ++path[level];
            continue;
        }

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

void mr_octree_periodic_wrap(mr_ocforest *forest, mr_index octree_idx, mr_float p[3]) {
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

mr_int mr_octree_locate_point(mr_ocforest *forest, mr_index octree_idx, const mr_float p[3]) {
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

        child_idx = (mr_int)(wp[0] > node->x) | (mr_int)(wp[1] > node->y) << 1 | (mr_int)(wp[2] > node->z) << 2;
        idx = node->first_child + child_idx;
    }

    if (mr_norm_inf(wp[0] - node->x, wp[1] - node->y, wp[2] - node->z) > node->dim / 2.0f) {
        return MR_INVALID_INDEX;
    }

    return idx;
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

typedef struct refinement_info {
    const mr_uint expected_level;

    bool did_refine;
    bool did_update_level;
} refinement_info;

static mr_int descend_path(
    mr_ocforest *forest,
    mr_int node_idx,
    mr_int path[MR_OCTREE_MAX_LEVEL],
    mr_uint len,
    refinement_info *info,
    mr_uint max_level
) {
    mr_int curr_idx = node_idx;
    mr_octree_node *node = NULL;

    while (len > 0) {
        node = mr_ocforest_get_node(forest, curr_idx);
        if (info && node->flags & MR_OCTREE_NODE_FLAG_LEAF && node->level < max_level) {
            insert_children(forest, curr_idx, NULL);
            node = mr_ocforest_get_node(forest, curr_idx);
            info->did_refine = true;
        }

        if (node->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            break;
        }

        curr_idx = node->first_child + path[--len];
    }

    return curr_idx;
}

static mr_int mr_octree_find_face_neighbor_with_refine(
    mr_ocforest *forest,
    mr_int idx,
    mr_direction dir,
    refinement_info *info
) {
    assert(forest);
    assert(idx != MR_INVALID_INDEX);

    bool cache_might_change = true;
    mr_int nidx = get_cached_neighbor(forest, idx, dir, &cache_might_change);
    if (nidx != MR_INVALID_INDEX && (!cache_might_change || !info)) {
        return nidx;
    }

    mr_uint len = 0;
    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };
    mr_uint max_level = info ? info->expected_level : 0;

    mr_int ancestor_idx = find_ancestor(forest, idx, &dir, MR_ADJACENCY_FACE, path, &len);
    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    if (info && len == 2) {
        ++max_level;
        info->did_update_level = true;
    }

    reflect_path(&dir, MR_ADJACENCY_FACE, path, len);
    nidx = descend_path(forest, ancestor_idx, path, len, info, max_level);
    cache_neighbor(forest, idx, nidx, dir);

    return nidx;
}

mr_int mr_octree_find_face_neighbor(mr_ocforest *forest, mr_int idx, mr_direction dir) {
    return mr_octree_find_face_neighbor_with_refine(forest, idx, dir, NULL);
}

static int common_face(
    mr_int local_idx,
    mr_direction dir[],
    mr_adjacency adj,
    mr_direction *face
) {
    if (adj < MR_ADJACENCY_EDGE) {
        return MR_FAILURE;
    }

    size_t cnt = 0;
    for (size_t i = 0; i < adj; ++i) {
        if (!is_node_adjacent(local_idx, &dir[i], MR_ADJACENCY_FACE)) {
            continue;
        }

        *face = dir[i];
        ++cnt;
    }

    if (cnt != 1) {
        return MR_FAILURE;
    }

    return MR_SUCCESS;
}

static mr_int mr_octree_find_edge_neighbor_with_refine(
    mr_ocforest *forest,
    mr_int idx,
    mr_direction dir[MR_ADJACENCY_EDGE],
    refinement_info *info
) {
    assert(forest);
    assert(dir);
    assert(idx != MR_INVALID_INDEX);

    mr_uint len = 0;
    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };
    mr_uint max_level = info ? info->expected_level : 0;

    mr_int ancestor_idx = find_ancestor(forest, idx, dir, MR_ADJACENCY_EDGE, path, &len);
    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_direction cface;
    mr_int neighbor_ancestor_idx = ancestor_idx;
    if (common_face(path[len - 1], dir, MR_ADJACENCY_EDGE, &cface) == MR_SUCCESS) {
        neighbor_ancestor_idx = mr_octree_find_face_neighbor_with_refine(forest, ancestor_idx, cface, info);
    }
    
    if (neighbor_ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    if (info && len == 2 && neighbor_ancestor_idx == ancestor_idx) {
        ++max_level;
        info->did_update_level = true;
    }

    reflect_path(dir, MR_ADJACENCY_EDGE, path, len);
    mr_int neighbor_idx = descend_path(forest, neighbor_ancestor_idx, path, len, info, max_level);
    if (neighbor_idx == neighbor_ancestor_idx) {
        return MR_INVALID_INDEX;
    }

    return neighbor_idx;
}

mr_int mr_octree_find_edge_neighbor(mr_ocforest *forest, mr_int idx, mr_direction dir[MR_ADJACENCY_EDGE]) {
    return mr_octree_find_edge_neighbor_with_refine(forest, idx, dir, NULL);
}

static int common_edge(
    mr_int local_idx,
    mr_direction dir[MR_ADJACENCY_VERTEX],
    mr_direction edge[MR_ADJACENCY_EDGE]
) {
    size_t cnt = 0;
    for (size_t i = 0; i < MR_ADJACENCY_VERTEX; ++i) {
        mr_direction curr_edge[] = { dir[i != 0 ? i - 1 : 2], dir[(i + 1) % 3] };
        if (!is_node_adjacent(local_idx, curr_edge, MR_ADJACENCY_EDGE)) {
            continue;
        }

        edge[0] = curr_edge[0];
        edge[1] = curr_edge[1];

        ++cnt;
    }

    if (cnt != 1) {
        return MR_FAILURE;
    }

    return MR_SUCCESS;
}

static mr_int mr_octree_find_vertex_neighbor_with_refine(
    mr_ocforest *forest,
    mr_int idx,
    mr_direction dir[MR_ADJACENCY_VERTEX],
    refinement_info *info
) {
    assert(forest);
    assert(dir);
    assert(idx != MR_INVALID_INDEX);

    mr_uint len = 0;
    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };
    mr_uint max_level = info ? info->expected_level : 0;

    mr_int ancestor_idx = find_ancestor(forest, idx, dir, MR_ADJACENCY_VERTEX, path, &len);
    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_direction common_elem[MR_ADJACENCY_EDGE];
    mr_int neighbor_ancestor_idx = ancestor_idx;
    if (common_edge(path[len - 1], dir, common_elem) == MR_SUCCESS) {
        neighbor_ancestor_idx = mr_octree_find_edge_neighbor_with_refine(forest, ancestor_idx, common_elem, info);
    } else if (common_face(path[len - 1], dir, MR_ADJACENCY_VERTEX, common_elem) == MR_SUCCESS) {
        neighbor_ancestor_idx = mr_octree_find_face_neighbor_with_refine(forest, ancestor_idx, common_elem[0], info);
    }

    if (neighbor_ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    if (info && len == 2 && neighbor_ancestor_idx == ancestor_idx) {
        ++max_level;
        info->did_update_level = true;
    }

    reflect_path(dir, MR_ADJACENCY_VERTEX, path, len);
    mr_int neighbor_idx = descend_path(forest, neighbor_ancestor_idx, path, len, info, max_level);
    if (neighbor_idx == neighbor_ancestor_idx) {
        return MR_INVALID_INDEX;
    }

    return neighbor_idx;
}

mr_int mr_octree_find_vertex_neighbor(mr_ocforest *forest, mr_int idx, mr_direction dir[MR_ADJACENCY_VERTEX]) {
    return mr_octree_find_vertex_neighbor_with_refine(forest, idx, dir, NULL);
}

static int activate_leaf(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);

    mr_octree_cond_cb *cond_ud = userdata;
    if (cond_ud && cond_ud->fn && !cond_ud->fn(forest, node_idx, cond_ud->userdata)) {
        node->flags &= ~MR_OCTREE_NODE_FLAG_ACTIVE;
    } else {
        node->flags |= MR_OCTREE_NODE_FLAG_ACTIVE;
    }

    return MR_SUCCESS;
}

void mr_octree_activate(mr_ocforest *forest, mr_index octree_idx, mr_octree_cond_cb cond) {
    mr_octree_leaves_apply(forest, octree_idx, mr_octree_cond_cb_null(), mr_octree_apply_cb_create(activate_leaf, &cond), false);
}

void mr_octree_activate_all(mr_ocforest *forest, mr_index octree_idx) {
    mr_octree_activate(forest, octree_idx, mr_octree_cond_cb_null());
}

static int insert_children(mr_ocforest *forest, mr_int parent_idx, void *userdata) {
    mr_octree_cond_cb *cond_ud = userdata;
    if (cond_ud && cond_ud->fn && !cond_ud->fn(forest, parent_idx, cond_ud->userdata)) {
        return MR_SUCCESS;
    }

    mr_octree_node *parent = mr_ocforest_get_node(forest, parent_idx);
    if (parent->level >= MR_OCTREE_MAX_LEVEL) {
        return MR_SUCCESS;
    }

    mr_index first_child = mr_mem_pool_alloc_many(forest->nodes, MR_OCTREE_NB_CHILDREN);
    if (first_child == MR_INVALID_INDEX) {
        return MR_FAILURE;
    }

    parent = mr_ocforest_get_node(forest, parent_idx);
    parent->flags &= ~MR_OCTREE_NODE_FLAG_LEAF;
    parent->first_child = first_child;

    mr_float qdim = parent->dim / 4.0f;
    for (size_t i = 0; i < MR_OCTREE_NB_CHILDREN; ++i) {
        mr_octree_node *node = mr_ocforest_get_node(forest, first_child + i);

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
            .first_child = MR_INVALID_INDEX,
            .value = parent->value,
            .gradient = { parent->gradient[0], parent->gradient[1], parent->gradient[2] },
        };

        mr_octree_node_connection *conn = mr_ocforest_get_node_connection(forest, first_child + i);
        *conn = (mr_octree_node_connection) {
            .local_idx = i,
            .external_neighbors = { MR_INVALID_INDEX, MR_INVALID_INDEX, MR_INVALID_INDEX },
        };
    }

    return MR_SUCCESS;
}

void mr_octree_refine(
    mr_ocforest *forest,
    mr_index octree_idx,
    mr_octree_cond_cb filter,
    mr_octree_cond_cb cond,
    bool recursive
) {
    mr_octree_leaves_apply(forest, octree_idx, filter, mr_octree_apply_cb_create(insert_children, &cond), recursive);
}

void mr_octree_refine_all(mr_ocforest *forest, mr_index octree_idx) {
    mr_octree_refine(forest, octree_idx, mr_octree_cond_cb_null(), mr_octree_cond_cb_null(), false);
}

typedef struct balance_userdata {
    mr_uint level;
    bool needs_repeat;
} balance_userdata;

static bool balance_filter(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    balance_userdata *balance_ud = userdata;
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);

    return node->level <= balance_ud->level;
}

static int balance_level(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    balance_userdata *balance_ud = userdata;

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    if (node->flags & (MR_OCTREE_NODE_FLAG_UNBALANCED | MR_OCTREE_NODE_FLAG_LEAF) && node->level == balance_ud->level) {
        mr_uint refine_level = node->level - 1;
        refinement_info info = { refine_level, false, false };

        for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
            mr_octree_find_face_neighbor_with_refine(forest, node_idx, dir, &info);
        }

        for (mr_direction dir_0 = MR_DIRECTION_MI_X; dir_0 <= MR_DIRECTION_PL_Y; ++dir_0) {
            for (mr_direction dir_1 = dir_0 - mr_direction_get_sign(dir_0) + 2; dir_1 <= MR_DIRECTION_PL_Z; ++dir_1) {
                mr_direction edge[] = { dir_0, dir_1 };
                mr_octree_find_edge_neighbor_with_refine(forest, node_idx, edge, &info);
            }
        }

// Vertex balancing isn't actually needed for current algorithms, but I'll leave it for now
#ifdef BALANCE_VERTEX
#define NB_VERTICES 8
        for (mr_uint i = 0; i < NB_VERTICES; ++i) {
            mr_sign signs[] = {
                i & 1,
                i >> 1 & 1,
                i >> 2 & 1,
            };

            mr_direction vertex[] = {
                mr_direction_create(MR_AXIS_X, signs[MR_AXIS_X]),
                mr_direction_create(MR_AXIS_Y, signs[MR_AXIS_Y]),
                mr_direction_create(MR_AXIS_Z, signs[MR_AXIS_Z]),
            };

            mr_octree_find_vertex_neighbor_with_refine(forest, node_idx, vertex, &info);
        }
#endif // BALANCE_VERTEX

        if (info.did_refine && info.did_update_level) {
            balance_ud->needs_repeat = true;
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
        balance_userdata userdata = {
            .level = i,
            .needs_repeat = false,
        };

        do {
            userdata.needs_repeat = false;

            mr_octree_leaves_apply(
                forest,
                octree_idx,
                mr_octree_cond_cb_create(balance_filter, &userdata),
                mr_octree_apply_cb_create(balance_level, &userdata),
                false
            );
        } while (userdata.needs_repeat);
    }
}