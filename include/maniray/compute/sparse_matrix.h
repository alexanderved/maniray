#ifndef _MR_SPARSE_MATRIX_H
#define _MR_SPARSE_MATRIX_H

#include <stdbool.h>

#include "maniray/utils/types.h"

typedef struct mr_sparse_row {
    size_t cap;
    size_t len;

    mr_float *values;
    size_t *cols;
} mr_sparse_row;

mr_sparse_row *mr_sparse_row_create();
void mr_sparse_row_destroy(mr_sparse_row *row);

mr_float mr_sparse_row_get(mr_sparse_row *row, size_t col);
void mr_sparse_row_set(mr_sparse_row *row, size_t col, mr_float val);

void mr_sparse_row_clear(mr_sparse_row *row);

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

mr_sparse_matrix_builder *mr_sparse_matrix_builder_copy(mr_sparse_matrix_builder *other);

void mr_sparse_matrix_builder_add_row(mr_sparse_matrix_builder *builder, mr_sparse_row *row);

#endif // _MR_SPARSE_MATRIX_H