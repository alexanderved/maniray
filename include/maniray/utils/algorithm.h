#ifndef _MR_ALGORITHM_H
#define _MR_ALGORITHM_H

#include <stdlib.h>

#include "maniray/utils/misc.h"
#include "maniray/utils/types.h"

// mr_index mr_binary

mr_index mr_generic_binary_search(
    const void *arr,
    size_t len,
    const void *key,
    size_t size,
    int (*cmp)(const void *, const void *)
);

void mr_generic_array_insert(
    void *arr,
    size_t end,
    const void *value,
    size_t size,
    mr_index pos
);

#define MR_GENERIC_CMP_NAME(type) mr_cmp_ ## type
#define MR_GENERIC_CMP_GEN(type) \
    static int MR_GENERIC_CMP_NAME(type)(const void *a, const void *b) { \
        return (*(type *)a > *(type *)b) - (*(type *)a < *(type *)b); \
    }

MR_GENERIC_CMP_GEN(mr_int)
MR_GENERIC_CMP_GEN(mr_uint)
MR_GENERIC_CMP_GEN(mr_float)
MR_GENERIC_CMP_GEN(mr_index)
MR_GENERIC_CMP_GEN(size_t)

#define MR_GENERIC_CMP_CHOOSE(key) \
    _Generic((key), \
        mr_int: MR_GENERIC_CMP_NAME(mr_int), \
        mr_uint: MR_GENERIC_CMP_NAME(mr_uint), \
        mr_float: MR_GENERIC_CMP_NAME(mr_float), \
        mr_index: MR_GENERIC_CMP_NAME(mr_index), \
        size_t: MR_GENERIC_CMP_NAME(size_t) \
    )

#define MR_BINARY_SEARCH(arr, len, key) \
    mr_generic_binary_search((arr), (len), &(key), sizeof((key)), MR_GENERIC_CMP_CHOOSE((key)))

#define MR_ARRAY_INSERT(arr, end, value, pos) \
    mr_generic_array_insert((arr), (end), &(value), sizeof((value)), (pos))

#endif // _MR_ALGORITHM_H