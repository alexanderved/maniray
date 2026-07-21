#ifndef _MR_OCTREE_H
#define _MR_OCTREE_H

#include "maniray/utils/types.h"
#include "maniray/utils/mem_pool.h"
#include "maniray/compute/manifold.h"

#define MR_OCTREE_MAX_LEVEL 6
#define MR_OCTREE_NB_CHILDREN 8

#define MR_OCTREE_NODE_FLAG_ACTIVE (1 << 0)
#define MR_OCTREE_NODE_FLAG_LEAF (1 << 1)

#define MR_OCTREE_FLAG_PERIODIC_X (1 << 0)
#define MR_OCTREE_FLAG_PERIODIC_Y (1 << 1)
#define MR_OCTREE_FLAG_PERIODIC_Z (1 << 2)

#define MR_OCTREE_NODE_FIELD 0

typedef enum mr_octree_direction {
    MR_OCTREE_DIRECTION_MI_X,
    MR_OCTREE_DIRECTION_PL_X,
    MR_OCTREE_DIRECTION_MI_Y,
    MR_OCTREE_DIRECTION_PL_Y,
    MR_OCTREE_DIRECTION_MI_Z,
    MR_OCTREE_DIRECTION_PL_Z,
} mr_octree_direction;

typedef struct mr_octree_node {
    mr_bitfield flags;
    mr_uint level;

    mr_uint chart_idx;
    mr_float x, y, z;
    mr_float dim;

    mr_int root;
    mr_int parent;
    mr_int first_child;

    mr_float value;
    mr_float gradient[3];
} mr_octree_node;

typedef struct mr_octree_root {
    mr_bitfield flags;
    mr_int node_idx;
} mr_octree_root;

typedef struct mr_octree_root_desc {
    mr_bitfield flags;

    mr_uint chart_idx;
    mr_float x, y, z;
    mr_float dim;
} mr_octree_root_desc;

typedef struct mr_ocforest {
    mr_manifold *manifold;

    size_t nb_roots;
    mr_octree_root *roots;

    mr_mem_pool *nodes;
} mr_ocforest;

typedef int (*mr_octree_apply_fn)(mr_ocforest *, mr_int, void *);
typedef struct mr_octree_apply_cb {
    mr_octree_apply_fn fn;
    void *userdata;
} mr_octree_apply_cb;

static inline mr_octree_apply_cb mr_make_octree_apply_cb(mr_octree_apply_fn fn, void *userdata) {
    return (mr_octree_apply_cb) { fn, userdata };
}

static inline mr_octree_apply_cb mr_octree_apply_cb_null() {
    return (mr_octree_apply_cb) { NULL, NULL };
}

typedef bool (*mr_octree_cond_fn)(mr_ocforest *, mr_int, void *);
typedef struct mr_octree_cond_cb {
    mr_octree_cond_fn fn;
    void *userdata;
} mr_octree_cond_cb;

static inline mr_octree_cond_cb mr_make_octree_cond_cb(mr_octree_cond_fn fn, void *userdata) {
    return (mr_octree_cond_cb) { fn, userdata };
}

static inline mr_octree_cond_cb mr_octree_cond_cb_null() {
    return (mr_octree_cond_cb) { NULL, NULL };
}

mr_ocforest *mr_ocforest_create(mr_manifold *manifold, mr_octree_root_desc roots[], size_t nb_roots);
void mr_ocforest_destroy(mr_ocforest *forest);

size_t mr_ocforest_size(mr_ocforest *forest);
mr_octree_node *mr_ocforest_get_node(mr_ocforest *forest, mr_int idx);

void mr_octree_apply(
    mr_ocforest *forest,
    mr_index octree_idx,
    mr_octree_cond_cb filter,
    mr_octree_apply_cb apply,
    bool recursive
);
void mr_octree_leaves_apply(
    mr_ocforest *forest,
    mr_index octree_idx,
    mr_octree_cond_cb filter,
    mr_octree_apply_cb apply,
    bool recursive
);

mr_int mr_octree_find_leaf(mr_ocforest *forest, mr_index octree_idx, mr_octree_cond_cb cond);
mr_int mr_octree_find_face_neighbor(mr_ocforest *forest, mr_int idx, mr_octree_direction dir);

void mr_octree_activate(mr_ocforest *forest, mr_index octree_idx, mr_octree_cond_cb cond);
void mr_octree_activate_all(mr_ocforest *forest, mr_index octree_idx);

void mr_octree_refine(
    mr_ocforest *forest,
    mr_index octree_idx,
    mr_octree_cond_cb filter,
    mr_octree_cond_cb cond,
    bool recursive
);
void mr_octree_refine_all(mr_ocforest *forest, mr_index octree_idx);

void mr_octree_balance(mr_ocforest *forest, mr_index octree_idx);

#endif // _MR_OCTREE_H