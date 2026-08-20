#ifndef _MR_OCTREE_H
#define _MR_OCTREE_H

#include "maniray/utils/mem_pool.h"
#include "maniray/compute/manifold.h"
#include "maniray/compute/navigation.h"

#define MR_OCTREE_NB_CHILDREN 8
#define MR_OCTREE_MAX_LEVEL 6

#define MR_OCTREE_NODE_BLOCK_DIM 4
#define MR_OCTREE_NB_CELLS_IN_BLOCK 64

#if MR_OCTREE_MAX_LEVEL > 6
#line 9
#error "The height of the octree is too large"
#endif

#define MR_OCTREE_NODE_FLAG_ACTIVE (1 << 0)
#define MR_OCTREE_NODE_FLAG_LEAF (1 << 1)
#define MR_OCTREE_NODE_FLAG_UNBALANCED (1 << 2)

#define MR_OCTREE_FLAG_PERIODIC_X (1 << MR_AXIS_X)
#define MR_OCTREE_FLAG_PERIODIC_Y (1 << MR_AXIS_Y)
#define MR_OCTREE_FLAG_PERIODIC_Z (1 << MR_AXIS_Z)

#define MR_OCTREE_NODE_NB_MAIN_FIELDS 2
#define MR_OCTREE_NODE_FIELD 0
#define MR_OCTREE_NODE_CONNECTION_FIELD 1

#define MR_OCTREE_CELL_NB_MAIN_FIELDS 1
#define MR_OCTREE_CELL_FIELD 0

typedef struct mr_octree_cell {
    mr_int parent;

    mr_uint chart_idx;
    mr_float x, y, z;
    mr_float dim;
} mr_octree_cell;

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

typedef struct mr_octree_node_connection {
    mr_int local_idx;
    mr_int external_neighbors[MR_NB_AXES];
} mr_octree_node_connection;

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
    mr_mem_pool *cells;
} mr_ocforest;

MR_DEFINE_CALLBACK(mr_octree_apply, int, mr_ocforest *forest, mr_int idx)
MR_DEFINE_CALLBACK(mr_octree_cond, bool, mr_ocforest *forest, mr_int idx)

mr_ocforest *mr_ocforest_create(
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots,
    size_t cell_extra_fields[],
    size_t nb_cell_extra_fields
);
void mr_ocforest_destroy(mr_ocforest *forest);

size_t mr_ocforest_nb_nodes_upper_bound(mr_ocforest *forest);
size_t mr_ocforest_nb_cells_upper_bound(mr_ocforest *forest);

size_t mr_ocforest_count_leaves(mr_ocforest *forest);
size_t mr_ocforest_count_cells(mr_ocforest *forest);

mr_octree_node *mr_ocforest_get_node(mr_ocforest *forest, mr_int idx);
mr_octree_node *mr_ocforest_get_node_array(mr_ocforest *forest);

mr_octree_node_connection *mr_ocforest_get_node_connection(mr_ocforest *forest, mr_int idx);
mr_octree_node_connection *mr_ocforest_get_node_connection_array(mr_ocforest *forest);

mr_octree_cell *mr_ocforest_get_cell(mr_ocforest *forest, mr_int idx);
mr_octree_cell *mr_ocforest_get_cell_array(mr_ocforest *forest);

void *mr_ocforest_get_cell_extra(mr_ocforest *forest, mr_int idx, mr_int field);
void *mr_ocforest_get_cell_extra_array(mr_ocforest *forest, mr_int field);

mr_int mr_ocforest_get_code(mr_ocforest *forest, mr_int cell_idx);
mr_int mr_ocforest_find_cell_with_code(mr_ocforest *forest, mr_int code);

int mr_octree_leaves_apply(mr_ocforest *forest, mr_index octree_idx, mr_octree_apply_cb apply);
int mr_octree_cells_apply(mr_ocforest *forest, mr_index octree_idx, mr_octree_apply_cb apply);

void mr_octree_periodic_wrap(mr_ocforest *forest, mr_index octree_idx, mr_float p[MR_NB_AXES]);
mr_int mr_octree_locate_point_in_leaf(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES]);
mr_int mr_octree_locate_point_in_cell(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES]);

mr_int mr_octree_find_face_neighbor_node(mr_ocforest *forest, mr_int idx, mr_direction dir);
mr_int mr_octree_find_face_neighbor_cell(mr_ocforest *forest, mr_int idx, mr_direction dir);

void mr_octree_refine(mr_ocforest *forest, mr_index octree_idx, mr_octree_cond_cb cond, bool recursive);
void mr_octree_refine_all(mr_ocforest *forest, mr_index octree_idx, size_t nb_repeats);

void mr_octree_balance(mr_ocforest *forest, mr_index octree_idx);


/*
 * TODO:
 * 1. Add userdata with dtor
 * 2. Add callbacks (call them mr_octree_ops) which should be called on some actions (activation, node initialization during refinement, etc.)
 * 3. Add functions to conviniently call this callbacks
 * 4. Add face neighbor caching
*/

#endif // _MR_OCTREE_H