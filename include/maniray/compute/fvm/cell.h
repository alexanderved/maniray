#ifndef _MR_FVM_CELL_H
#define _MR_FVM_CELL_H

#include "maniray/utils/types.h"
#include "maniray/compute/octree.h"

mr_float mr_cell_volume(mr_ocforest *forest, mr_int idx);
mr_float mr_cell_face_area(mr_ocforest *forest, mr_int idx, mr_direction face_dir);
mr_float mr_cell_neighbor_distance(mr_ocforest *forest, mr_int idx, mr_direction dir);

#endif // _MR_FVM_CELL_H