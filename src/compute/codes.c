#include "maniray/utils/algorithm.h"
#include "maniray/utils/xmalloc.h"
#include "maniray/compute/codes.h"

typedef struct code_map_userdata {
    mr_int curr_idx;
    mr_int *codes;
} code_map_userdata;

static int save_code(mr_ocforest *forest, mr_int node_idx, void *userdata) {
    code_map_userdata *code_map_ud = userdata;
    code_map_ud->codes[code_map_ud->curr_idx++] = mr_ocforest_get_code(forest, node_idx);

    return MR_SUCCESS;
}

mr_code_map *mr_code_map_create_from_ocforest(mr_ocforest *forest) {
    mr_code_map *map = xmalloc(sizeof(mr_code_map));

    size_t nb_leaves = mr_ocforest_count_leaves(forest);
    map->len = nb_leaves;
    map->codes = xmalloc(nb_leaves * sizeof(mr_int));

    code_map_userdata ud = { 0, map->codes };
    for (mr_index octree_idx = 0; octree_idx < forest->nb_roots; ++octree_idx) {
        mr_octree_leaves_apply(forest, octree_idx, mr_octree_cond_cb_null(), mr_octree_apply_cb_create(save_code, &ud), false);
    }

    return map;
}

void mr_code_map_destroy(mr_code_map *map) {
    if (!map) {
        return;
    }

    free(map->codes);
    free(map);
}

mr_int mr_code_map_get_code(mr_code_map *map, mr_int idx) {
    return map && idx < map->len ? map->codes[idx] : MR_INVALID_INDEX;
}

mr_int mr_code_map_get_index(mr_code_map *map, mr_int code) {
    return map ? MR_BINARY_SEARCH(map->codes, map->len, code) : MR_INVALID_INDEX;
}