#include <stddef.h>
#include <stdalign.h>
#include <string.h>

#include "maniray/utils/dyn_array.h"
#include "maniray/utils/xmalloc.h"

typedef struct array_header {
    alignas(alignof(max_align_t)) size_t type_size;
    size_t cap;
    size_t len;
} array_header;

void *mr_dyn_array_create_with_info(size_t type_size, size_t cap, size_t len, const void *data) {
    if (len > cap) {
        cap = len;
    }

    array_header *header = xmalloc(sizeof(array_header) + cap * type_size);
    header->type_size = type_size;
    header->cap = cap;
    header->len = len;

    if (len != 0 && data) {
        memcpy(header + 1, data, len);
    }

    return header + 1;
}

void mr_dyn_array_destroy(void *arr) {
    if (!arr) {
        return;
    }

    array_header *header = (array_header *)arr - 1;
    free(header);
}

void *mr_dyn_array_push_alloc(void *arr, void *elem) {
    array_header *header = (array_header *)arr - 1;
    if (header->len == header->cap) {
        if (header->cap == 0) {
            header->cap = 1;
        } else {
            header->cap *= 2;
        }

        header = xrealloc(header, sizeof(array_header) + header->cap * header->type_size);
    }

    char *arr_elem = (char *)header + sizeof(array_header) + header->len * header->type_size;
    memcpy(arr_elem, elem, header->type_size);
    ++header->len;

    return header + 1;
}

void mr_dyn_array_pop(void *arr) {
    array_header *header = (array_header *)arr - 1;
    --header->len;
}