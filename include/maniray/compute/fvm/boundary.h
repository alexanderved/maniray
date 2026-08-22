#ifndef _MR_FVM_BOUNDARY_H
#define _MR_FVM_BOUNDARY_H

#include <stdbool.h>

#include "maniray/compute/octree.h"

typedef struct mr_boundary_condition {
    int _currently_nothing;
} mr_boundary_condition;

bool mr_is_boundary_cell(mr_ocforest *forest, mr_int cell_idx, mr_direction dir);

#endif // _MR_FVM_BOUNDARY_H