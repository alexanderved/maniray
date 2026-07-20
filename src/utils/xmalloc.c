#include <stdlib.h>

#include "maniray/utils/xmalloc.h"

void *xmalloc(size_t size) {
    if (size == 0) {
        exit(1);
    }

    void *ptr = malloc(size);
    if (!ptr) {
        exit(1);
    }

    return ptr;
}

void *xcalloc(size_t n, size_t size) {
    if (size == 0) {
        exit(1);
    }

    void *ptr = calloc(n, size);
    if (!ptr) {
        exit(1);
    }

    return ptr;
}

void *xrealloc(void *ptr, size_t size) {
    if (size == 0) {
        free(ptr);
        exit(1);
    }

    void *new_ptr = realloc(ptr, size);
    if (!new_ptr) {
        free(ptr);
        exit(1);
    }

    return new_ptr;
}