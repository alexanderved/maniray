#include <string.h>
#include <math.h>
#include <stdio.h>

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

mr_octree_node *mr_ocforest_get_node(mr_ocforest *forest, mr_int idx) {
    return forest ? mr_mem_pool_ptr(forest->nodes, MR_OCTREE_NODE_FIELD, idx) : NULL;
}

mr_octree_node *mr_ocforest_get_node_array(mr_ocforest *forest) {
    return forest ? mr_mem_pool_array_ptr(forest->nodes, MR_OCTREE_NODE_FIELD) : NULL;
}

void *mr_ocforest_get_extra(mr_ocforest *forest, mr_int idx, mr_int field) {
    return forest ? mr_mem_pool_ptr(forest->nodes, MR_OCTREE_NODE_NB_MAIN_FIELDS + field, idx) : NULL;
}

void *mr_ocforest_get_extra_array(mr_ocforest *forest, mr_int field) {
    return forest ? mr_mem_pool_array_ptr(forest->nodes, MR_OCTREE_NODE_NB_MAIN_FIELDS + field) : NULL;
}

int mr_octree_leaves_apply(
    mr_ocforest *forest,
    mr_index octree_idx,
    mr_octree_cond_cb filter,
    mr_octree_apply_cb apply,
    bool recursive
) {
    if (!forest || (size_t)octree_idx >= forest->nb_roots) {
        return MR_FAILURE;
    }

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

static mr_int get_local_idx(mr_ocforest *forest, mr_int node_idx) {
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    mr_octree_node *parent = mr_ocforest_get_node(forest, node->parent);
    if (!parent) {
        return MR_INVALID_INDEX;
    }

    return node_idx - parent->first_child;
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
    size_t *nb_refine,
    mr_uint max_level
) {
    mr_int curr_idx = node_idx;
    mr_octree_node *node = NULL;

    while (len > 0) {
        node = mr_ocforest_get_node(forest, curr_idx);
        if (nb_refine && node->flags & MR_OCTREE_NODE_FLAG_LEAF && node->level < max_level) {
            insert_children(forest, curr_idx, NULL);
            node = mr_ocforest_get_node(forest, curr_idx);
            ++(*nb_refine);
        }

        if (node->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            break;
        }

        curr_idx = node->first_child + path[--len];
    }

    return curr_idx;
}

static mr_int mr_octree_find_face_neighbor_ext(mr_ocforest *forest, mr_int idx, mr_direction dir, size_t *nb_refine, mr_uint max_level) {
    if (!forest || idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_uint len = 0;
    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };

    mr_int ancestor_idx = find_ancestor(forest, idx, &dir, MR_ADJACENCY_FACE, path, &len);
    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    reflect_path(&dir, MR_ADJACENCY_FACE, path, len);

    return descend_path(forest, ancestor_idx, path, len, nb_refine, max_level);
}

mr_int mr_octree_find_face_neighbor(mr_ocforest *forest, mr_int idx, mr_direction dir) {
    return mr_octree_find_face_neighbor_ext(forest, idx, dir, NULL, 0);
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

static mr_int mr_octree_find_edge_neighbor_ext(
    mr_ocforest *forest,
    mr_int idx,
    mr_direction dir[MR_ADJACENCY_EDGE],
    size_t *nb_refine,
    mr_uint max_level
) {
    if (!forest || idx == MR_INVALID_INDEX || !dir) {
        return MR_INVALID_INDEX;
    }

    mr_uint len = 0;
    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };

    mr_int ancestor_idx = find_ancestor(forest, idx, dir, MR_ADJACENCY_EDGE, path, &len);
    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_direction cface;
    if (common_face(path[len - 1], dir, MR_ADJACENCY_EDGE, &cface) == MR_SUCCESS) {
        ancestor_idx = mr_octree_find_face_neighbor_ext(forest, ancestor_idx, cface, nb_refine, max_level);
    }
    
    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    reflect_path(dir, MR_ADJACENCY_EDGE, path, len);
    mr_int neighbor_idx = descend_path(forest, ancestor_idx, path, len, nb_refine, max_level);
    if (neighbor_idx == ancestor_idx) {
        return MR_INVALID_INDEX;
    }

    return neighbor_idx;
}

mr_int mr_octree_find_edge_neighbor(mr_ocforest *forest, mr_int idx, mr_direction dir[MR_ADJACENCY_EDGE]) {
    return mr_octree_find_edge_neighbor_ext(forest, idx, dir, NULL, 0);
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

static mr_int mr_octree_find_vertex_neighbor_ext(
    mr_ocforest *forest,
    mr_int idx,
    mr_direction dir[MR_ADJACENCY_VERTEX],
    size_t *nb_refine,
    mr_uint max_level
) {
    if (!forest || idx == MR_INVALID_INDEX || !dir) {
        return MR_INVALID_INDEX;
    }

    mr_uint len = 0;
    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };

    mr_int ancestor_idx = find_ancestor(forest, idx, dir, MR_ADJACENCY_VERTEX, path, &len);
    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_direction common_elem[MR_ADJACENCY_EDGE];
    if (common_edge(path[len - 1], dir, common_elem) == MR_SUCCESS) {
        ancestor_idx = mr_octree_find_edge_neighbor_ext(forest, ancestor_idx, common_elem, nb_refine, max_level);
    } else if (common_face(path[len - 1], dir, MR_ADJACENCY_VERTEX, common_elem) == MR_SUCCESS) {
        ancestor_idx = mr_octree_find_face_neighbor_ext(forest, ancestor_idx, common_elem[0], nb_refine, max_level);
    }

    if (ancestor_idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    reflect_path(dir, MR_ADJACENCY_VERTEX, path, len);
    mr_int neighbor_idx = descend_path(forest, ancestor_idx, path, len, nb_refine, max_level);
    if (neighbor_idx == ancestor_idx) {
        return MR_INVALID_INDEX;
    }

    return neighbor_idx;
}

mr_int mr_octree_find_vertex_neighbor(mr_ocforest *forest, mr_int idx, mr_direction dir[MR_ADJACENCY_VERTEX]) {
    return mr_octree_find_vertex_neighbor_ext(forest, idx, dir, false, 0);
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
            .flags = MR_OCTREE_NODE_FLAG_LEAF | (parent->flags & MR_OCTREE_NODE_FLAG_ACTIVE),
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
        };

        for (size_t i = 0; i < 3; ++i) {
            node->gradient[i] = parent->gradient[i];
        }
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

static bool is_mixed_node(mr_ocforest *forest, mr_octree_node *node) {
    if (node->parent == MR_INVALID_INDEX) {
        return false;
    }

    mr_octree_node *parent = mr_ocforest_get_node(forest, node->parent);
    for (mr_int i = 0; i < MR_OCTREE_NB_CHILDREN; ++i) {
        mr_octree_node *child = mr_ocforest_get_node(forest, parent->first_child + i);
        if (!(child->flags & MR_OCTREE_NODE_FLAG_LEAF)) {
            return true;
        }
    }

    return false;
}

static int balance_level(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    balance_userdata *balance_ud = userdata;

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    if (node->flags & MR_OCTREE_NODE_FLAG_LEAF && node->level == balance_ud->level) {
        bool is_mixed = is_mixed_node(forest, node);
        mr_uint refine_level;
        if (is_mixed) {
            refine_level = node->level;
        } else {
            refine_level = node->level - 1;
        }

        size_t nb_refine = 0;


        for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
            mr_octree_find_face_neighbor_ext(forest, node_idx, dir, &nb_refine, refine_level);
        }

        for (mr_direction dir_0 = MR_DIRECTION_MI_X; dir_0 <= MR_DIRECTION_PL_Y; ++dir_0) {
            for (mr_direction dir_1 = dir_0 - mr_direction_get_sign(dir_0) + 2; dir_1 <= MR_DIRECTION_PL_Z; ++dir_1) {
                mr_direction edge[] = { dir_0, dir_1 };
                mr_octree_find_edge_neighbor_ext(forest, node_idx, edge, &nb_refine, refine_level);
            }
        }

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

            mr_octree_find_vertex_neighbor_ext(forest, node_idx, vertex, &nb_refine, refine_level);
        }

        if (is_mixed && nb_refine != 0) {
            balance_ud->needs_repeat = true;
        }
    }

    return MR_SUCCESS;
}

void mr_octree_balance(mr_ocforest *forest, mr_index octree_idx) {
    if (!forest || (size_t)octree_idx >= forest->nb_roots) {
        return;
    }

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