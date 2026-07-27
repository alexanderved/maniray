#include "string.h"

#include "maniray/compute/sparse_matrix.h"
#include "maniray/utils/xmalloc.h"

mr_sparse_matrix_builder *mr_sparse_matrix_builder_create(size_t dim) {
    mr_sparse_matrix_builder *builder = xmalloc(sizeof(mr_sparse_matrix_builder));

    builder->dim = dim;
    builder->nb_cols = 0;
    builder->nb_rows = 0;

    builder->values = NULL;
    builder->cols = NULL;

    builder->rows = xmalloc((builder->dim + 1) * sizeof(size_t));
    builder->rows[0] = 0;

    return builder;
}

void mr_sparse_matrix_builder_destroy(mr_sparse_matrix_builder *builder) {
    free(builder->rows);
    free(builder->cols);
    free(builder->values);

    free(builder);
}

void mr_sparse_matrix_builder_add_row(mr_sparse_matrix_builder *builder, mr_float *values, size_t *cols, size_t len) {
    if (!builder || !values || !cols || builder->nb_rows == builder->dim) {
        return;
    }

    size_t old_nb_cols = builder->nb_cols;
    builder->nb_cols += len;
    builder->rows[++builder->nb_rows] = builder->nb_cols;

    builder->values = xrealloc(builder->values, builder->nb_cols * sizeof(mr_float));
    memcpy(builder->values + old_nb_cols, values, len * sizeof(mr_float));

    builder->cols = xrealloc(builder->cols, builder->nb_cols * sizeof(size_t));
    memcpy(builder->cols + old_nb_cols, cols, len * sizeof(size_t));
}
