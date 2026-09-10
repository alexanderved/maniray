#ifndef _MR_GEOMETRY_H
#define _MR_GEOMETRY_H

#include "maniray/utils/types.h"
#include "maniray/compute/octree.h"

void mr_cell_face_center(mr_ocforest *forest, mr_int idx, mr_direction face_dir, mr_float face_center[MR_NB_AXES]);

mr_float mr_cell_center_metric(mr_ocforest *forest, mr_int idx, size_t i, size_t j);
mr_float mr_cell_face_center_metric(mr_ocforest *forest, mr_int idx, mr_direction face_dir, size_t i, size_t j);

mr_float mr_cell_center_inv_metric(mr_ocforest *forest, mr_int idx, size_t i, size_t j);
mr_float mr_cell_face_center_inv_metric(mr_ocforest *forest, mr_int idx, mr_direction face_dir, size_t i, size_t j);

mr_float mr_cell_volume(mr_ocforest *forest, mr_int idx);
mr_float mr_cell_face_area(mr_ocforest *forest, mr_int idx, mr_direction face_dir);

#endif // _MR_GEOMETRY_H