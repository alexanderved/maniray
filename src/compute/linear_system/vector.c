#include <stdlib.h>

#include "maniray/compute/linear_system/vector.h"
#include "maniray/utils/xmalloc.h"

#include <lis_vector.h>

mr_vector *mr_vector_create(mr_float64 *arr, mr_int len) {
    mr_vector *vec = xmalloc(sizeof(mr_vector));

    vec->len = len;
    vec->data = arr;

    lis_vector_create(LIS_COMM_WORLD, &vec->inner);
    lis_vector_set_size(vec->inner, 0, vec->len);

    lis_vector_set(vec->inner, vec->data);
    vec->inner->status = LIS_VECTOR_ASSEMBLED;
    vec->inner->is_destroy = LIS_FALSE;

    return vec;
}

void mr_vector_destroy(mr_vector *vec) {
    if (!vec) {
        return;
    }

    lis_vector_destroy(vec->inner);
    free(vec->data);
    free(vec);
}