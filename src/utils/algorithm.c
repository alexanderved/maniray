#include <assert.h>
#include <string.h>

#include "maniray/utils/algorithm.h"

static inline const void *get_elem_ptr(const void *arr, size_t size, mr_index i) {
    return (const char *)arr + i * size;
}

mr_index mr_generic_binary_search(
    const void *arr,
    size_t len,
    const void *key,
    size_t size,
    int (*cmp)(const void *, const void *)
) {
    if (!arr || !key || !cmp) {
        return MR_INVALID_INDEX;
    }

    mr_index low = 0;
    mr_index high = len - 1;

    while (low <= high) {
        mr_index mid = low + (high - low) / 2;
        if (cmp(get_elem_ptr(arr, size, mid), key) < 0) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }

    return low;
}

void mr_generic_array_insert(
    void *arr,
    size_t end,
    const void *value,
    size_t size,
    mr_index pos
) {
    assert(arr);
    assert(value);

    mr_index offset = pos * size;
    char *ptr = (char *)arr + offset;

    if ((size_t)pos != end) {
        memmove(ptr + size, ptr, (end - pos) * size);
    }

    memcpy(ptr, value, size);
}