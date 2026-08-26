#ifndef _MR_FVM_BOUNDARY_H
#define _MR_FVM_BOUNDARY_H

#include <stdbool.h>

#include "maniray/compute/octree.h"

bool mr_is_boundary_cell(mr_ocforest *forest, mr_int cell_idx, mr_direction dir);

typedef enum mr_boundary_condition_type {
    MR_BC_NONE,
    MR_BC_DIRICHLET,
    MR_BC_NEUMANN,
} mr_boundary_condition_type;

typedef struct mr_boundary_condition mr_boundary_condition;
typedef void (*mr_boundary_condition_fn)(mr_boundary_condition *bc, mr_int cell_idx, mr_direction dir, void *out);

struct mr_boundary_condition {
    mr_ocforest *forest;
    mr_boundary_condition_type *types;
    mr_boundary_condition_fn *fns;
};

mr_boundary_condition *mr_boundary_condition_create(
    mr_ocforest *forest,
    mr_boundary_condition_type types[][MR_NB_DIRECTIONS],
    mr_boundary_condition_fn fns[]
);
void mr_boundary_condition_destroy(mr_boundary_condition *bc);

mr_boundary_condition_type mr_boundary_condition_get_type(mr_boundary_condition *bc, mr_int cell_idx, mr_direction dir);
void mr_boundary_condition_get_value(mr_boundary_condition *bc, mr_int cell_idx, mr_direction dir, void *out);

#endif // _MR_FVM_BOUNDARY_H