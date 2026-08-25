#ifndef _MR_LIN_SYS_VECTOR_H
#define _MR_LIN_SYS_VECTOR_H

#include <lis.h>

#include "maniray/utils/types.h"

typedef struct mr_vector {
    mr_int len;
    mr_float64 *data;

    LIS_VECTOR inner;
} mr_vector;

mr_vector *mr_vector_create(mr_float64 *arr, mr_int len);
void mr_vector_destroy(mr_vector *vec);

#endif // _MR_LIN_SYS_VECTOR_H