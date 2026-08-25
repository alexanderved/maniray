#ifndef _MR_LIN_SYS_MATRIX_H
#define _MR_LIN_SYS_MATRIX_H

#include <stdbool.h>

#include <lis.h>

#include "maniray/utils/types.h"

typedef struct mr_sparse_row {
    size_t cap;
    size_t len;

    mr_float64 *values;
    mr_int *cols;
} mr_sparse_row;

mr_sparse_row *mr_sparse_row_create();
void mr_sparse_row_destroy(mr_sparse_row *row);

mr_float64 mr_sparse_row_get(mr_sparse_row *row, mr_int col);
void mr_sparse_row_set(mr_sparse_row *row, mr_int col, mr_float64 val);

void mr_sparse_row_clear(mr_sparse_row *row);

typedef struct mr_sparse_matrix_builder {
    mr_int dim;
    mr_int nb_cols;
    mr_int nb_rows;

    mr_float64 *values;
    mr_int *cols;
    mr_int *rows;
} mr_sparse_matrix_builder;

mr_sparse_matrix_builder *mr_sparse_matrix_builder_create(mr_int dim);
void mr_sparse_matrix_builder_destroy(mr_sparse_matrix_builder *builder);

void mr_sparse_matrix_builder_add_row(mr_sparse_matrix_builder *builder, mr_sparse_row *row);

typedef struct mr_sparse_matrix {
    mr_int dim;

    mr_float64 *values;
    mr_int *cols;
    mr_int *rows;

    LIS_MATRIX inner;
} mr_sparse_matrix;

mr_sparse_matrix *mr_sparse_matrix_build(mr_sparse_matrix_builder *builder);
void mr_sparse_matrix_destroy(mr_sparse_matrix *mat);

#endif // _MR_LIN_SYS_MATRIX_H