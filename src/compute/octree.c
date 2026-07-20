#include "maniray/compute/octree.h"
#include "maniray/utils/xmalloc.h"
#include "maniray/utils/misc.h"

mr_ocforest *mr_ocforest_create(mr_manifold *manifold, mr_octree_root_desc roots[], size_t nb_roots) {
    mr_ocforest *forest = xmalloc(sizeof(mr_ocforest));

    forest->manifold = manifold;
    forest->nb_roots = nb_roots;
    forest->roots = xmalloc(nb_roots * sizeof(mr_octree_root));
    forest->nodes = mr_mem_pool_create((size_t[]) { sizeof(mr_octree_node) }, 1);

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

        mr_octree_node *node = mr_mem_pool_ptr(forest->nodes, 0, node_idx + i);
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

void mr_octree_apply(
    mr_ocforest *forest,
    mr_index octree_idx,
    mr_octree_cond_cb filter,
    mr_octree_apply_cb apply
) {
    if (!forest) {
        return;
    }

    mr_int root_idx = forest->roots[octree_idx].node_idx;
    mr_octree_node *root = mr_mem_pool_ptr(forest->nodes, 0, root_idx);
    if (root->first_child == MR_INVALID_INDEX) {
        apply.fn(forest, root_idx, apply.userdata);

        return;
    }

    mr_uint parent_idx = root_idx;
    mr_octree_node *parent = NULL;

    mr_uint child_idx = 0;
    mr_octree_node *child = NULL;

    mr_uint path[MR_OCTREE_MAX_LEVEL] = { 0 };
    mr_int level = 0;
    while (true) {
        parent = mr_mem_pool_ptr(forest->nodes, 0, parent_idx);
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

        child = mr_mem_pool_ptr(forest->nodes, 0, child_idx);

        if (child->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            apply.fn(forest, child_idx, apply.userdata);
            ++path[level];
        } else {
            parent_idx = child_idx;
            path[++level] = 0;
        }
    }
}

static int activate_leaf(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    mr_octree_cond_cb *cond_ud = userdata;
    if (cond_ud && cond_ud->fn && !cond_ud->fn(forest, node_idx, cond_ud->userdata)) {
        return MR_SUCCESS;
    }

    mr_octree_node *node = mr_mem_pool_ptr(forest->nodes, 0, node_idx);
    node->flags |= MR_OCTREE_NODE_FLAG_ACTIVE;

    return MR_SUCCESS;
}

void mr_octree_activate_all(mr_ocforest *forest, mr_index octree_idx) {
    mr_octree_apply(forest, octree_idx, mr_octree_cond_cb_null(), mr_make_octree_apply_cb(activate_leaf, NULL));
}

static int insert_children(mr_ocforest *forest, mr_int parent_idx, void *userdata) {
    mr_octree_cond_cb *cond_ud = userdata;
    if (cond_ud && cond_ud->fn && !cond_ud->fn(forest, parent_idx, cond_ud->userdata)) {
        return MR_SUCCESS;
    }

    mr_octree_node *parent = mr_mem_pool_ptr(forest->nodes, 0, parent_idx);
    if (parent->level >= MR_OCTREE_MAX_LEVEL) {
        return MR_SUCCESS;
    }

    mr_index first_child = mr_mem_pool_alloc_many(forest->nodes, MR_OCTREE_NB_CHILDREN);
    if (first_child == MR_INVALID_INDEX) {
        return MR_FAILURE;
    }

    parent = mr_mem_pool_ptr(forest->nodes, 0, parent_idx);
    parent->flags &= ~MR_OCTREE_NODE_FLAG_LEAF;
    parent->first_child = first_child;

    mr_float qdim = parent->dim / 4.0f;
    for (size_t i = 0; i < MR_OCTREE_NB_CHILDREN; ++i) {
        mr_octree_node *node = mr_mem_pool_ptr(forest->nodes, 0, first_child + i);

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
            .value = parent->value / MR_OCTREE_NB_CHILDREN,
        };

        for (size_t i = 0; i < 3; ++i) {
            node->gradient[i] = parent->gradient[i];
        }
    }

    return MR_SUCCESS;
}

void mr_octree_refine_all(mr_ocforest *forest, mr_index octree_idx) {
    mr_octree_refine(forest, octree_idx, mr_octree_cond_cb_null(), mr_octree_cond_cb_null());
}

void mr_octree_refine(
    mr_ocforest *forest,
    mr_index octree_idx,
    mr_octree_cond_cb filter,
    mr_octree_cond_cb cond
) {
    mr_octree_apply(forest, octree_idx, filter, mr_make_octree_apply_cb(insert_children, &cond));
}