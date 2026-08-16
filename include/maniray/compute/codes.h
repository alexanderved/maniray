#ifndef _MR_CODES_H
#define _MR_CODES_H

#include <stdlib.h>

#include "maniray/utils/types.h"
#include "maniray/compute/octree.h"

typedef struct mr_code_map {
    size_t len;
    mr_int *codes;
} mr_code_map;

mr_code_map *mr_code_map_create_from_ocforest(mr_ocforest *forest);
void mr_code_map_destroy(mr_code_map *map);

mr_int mr_code_map_get_code(mr_code_map *map, mr_int idx);
mr_int mr_code_map_get_index(mr_code_map *map, mr_int code);

#endif // _MR_CODES_H