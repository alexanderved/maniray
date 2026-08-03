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

mr_int mr_octree_locate_point(mr_ocforest *forest, mr_index octree_idx, mr_float p[3]) {
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

static mr_int mr_octree_find_face_neighbor_ext(mr_ocforest *forest, mr_int idx, mr_octree_direction dir, bool do_refine, mr_uint max_level) {
    if (!forest || idx == MR_INVALID_INDEX) {
        return MR_INVALID_INDEX;
    }

    mr_int dir_shift = dir / 2;
    mr_int dir_state = dir % 2;

    mr_octree_node *node = NULL;
    mr_octree_node *parent = NULL;

    mr_int curr_idx = idx;
    mr_int local_idx = 0;

    mr_uint len = 0;
    mr_int path[MR_OCTREE_MAX_LEVEL] = { 0 };
    do {
        node = mr_ocforest_get_node(forest, curr_idx);
        parent = mr_ocforest_get_node(forest, node->parent);
        if (!parent) {
            mr_octree_root *root = &forest->roots[node->root];
            if (((root->flags & MR_OCTREE_FLAG_PERIODIC_X) && (dir == MR_OCTREE_DIRECTION_MI_X || dir == MR_OCTREE_DIRECTION_PL_X))
                || ((root->flags & MR_OCTREE_FLAG_PERIODIC_Y) && (dir == MR_OCTREE_DIRECTION_MI_Y || dir == MR_OCTREE_DIRECTION_PL_Y))
                || ((root->flags & MR_OCTREE_FLAG_PERIODIC_Z) && (dir == MR_OCTREE_DIRECTION_MI_Z || dir == MR_OCTREE_DIRECTION_PL_Z)))
            {
                break;
            }

            return MR_INVALID_INDEX;
        }

        local_idx = curr_idx - parent->first_child;

        path[len++] = local_idx;
        curr_idx = node->parent;
    } while ((local_idx >> dir_shift & 1) == dir_state);

    do {
        parent = mr_ocforest_get_node(forest, curr_idx);

        local_idx = path[--len] ^ (1 << dir_shift);
        curr_idx = parent->first_child + local_idx;
        node = mr_ocforest_get_node(forest, curr_idx);

        if (do_refine && node->flags & MR_OCTREE_NODE_FLAG_LEAF && node->level < max_level) {
            insert_children(forest, curr_idx, NULL);
            node = mr_ocforest_get_node(forest, curr_idx);
        }
    } while (len > 0 && !(node->flags & MR_OCTREE_NODE_FLAG_LEAF));

    return curr_idx;
}

mr_int mr_octree_find_face_neighbor(mr_ocforest *forest, mr_int idx, mr_octree_direction dir) {
    return mr_octree_find_face_neighbor_ext(forest, idx, dir, false, 0);
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
    mr_octree_leaves_apply(forest, octree_idx, mr_octree_cond_cb_null(), mr_make_octree_apply_cb(activate_leaf, &cond), false);
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
    mr_octree_leaves_apply(forest, octree_idx, filter, mr_make_octree_apply_cb(insert_children, &cond), recursive);
}

void mr_octree_refine_all(mr_ocforest *forest, mr_index octree_idx) {
    mr_octree_refine(forest, octree_idx, mr_octree_cond_cb_null(), mr_octree_cond_cb_null(), false);
}

static bool balance_filter(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    mr_uint level = *(mr_int *)userdata;
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);

    return node->level <= level;
}

static int balance_level(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    mr_uint level = *(mr_int *)userdata;

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    if (node->flags & MR_OCTREE_NODE_FLAG_LEAF && node->level == level) {
        mr_uint node_level = node->level;
        for (mr_octree_direction dir = MR_OCTREE_DIRECTION_MI_X; dir <= MR_OCTREE_DIRECTION_PL_Z; ++dir) {
            mr_octree_find_face_neighbor_ext(forest, node_idx, dir, true, node_level - 1);
        }
    }

    return MR_SUCCESS;
}

void mr_octree_balance(mr_ocforest *forest, mr_index octree_idx) {
    if (!forest || (size_t)octree_idx >= forest->nb_roots) {
        return;
    }

    for (mr_uint i = MR_OCTREE_MAX_LEVEL; i > 0; --i) {
        mr_octree_leaves_apply(
            forest,
            octree_idx,
            mr_make_octree_cond_cb(balance_filter, &i),
            mr_make_octree_apply_cb(balance_level, &i),
            false
        );
    }
}