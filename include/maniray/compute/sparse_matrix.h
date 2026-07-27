#ifndef _MR_SPARSE_MATRIX_H
#define _MR_SPARSE_MATRIX_H

#include <stdbool.h>

#include "maniray/utils/types.h"

typedef struct mr_sparse_matrix_builder {
    size_t dim;
    size_t nb_cols;
    size_t nb_rows;

    mr_float *values;
    size_t *cols;
    size_t *rows;
} mr_sparse_matrix_builder;

mr_sparse_matrix_builder *mr_sparse_matrix_builder_create(size_t dim);
void mr_sparse_matrix_builder_destroy(mr_sparse_matrix_builder *builder);

void mr_sparse_matrix_builder_add_row(mr_sparse_matrix_builder *builder, mr_float *values, size_t *cols, size_t len);

#endif // _MR_SPARSE_MATRIX_H