#include <string.h>

#include "maniray/compute/sparse_matrix.h"
#include "maniray/utils/xmalloc.h"
#include "maniray/utils/algorithm.h"

mr_sparse_row *mr_sparse_row_create() {
    mr_sparse_row *row = xmalloc(sizeof(mr_sparse_row));

    row->cap = 0;
    row->len = 0;

    row->cols = NULL;
    row->values = NULL;

    return row;
}

void mr_sparse_row_destroy(mr_sparse_row *row) {
    if (!row) {
        return;
    }

    free(row->values);
    free(row->cols);

    free(row);
}

static void realloc_row(mr_sparse_row *row) {
    if (row->cap == 0) {
        ++row->cap;
        row->cols = xmalloc(sizeof(size_t));
        row->values = xmalloc(sizeof(mr_float));

        return;
    }

    row->cap *= 2;
    row->cols = xrealloc(row->cols, row->cap * sizeof(size_t));
    row->values = xrealloc(row->values, row->cap * sizeof(mr_float));
}

static void row_insert(mr_sparse_row *row, size_t col, mr_float val, mr_index pos) {
    if (row->len == row->cap) {
        realloc_row(row);
    }

    MR_ARRAY_INSERT(row->cols, row->len, col, pos);
    MR_ARRAY_INSERT(row->values, row->len, val, pos);

    ++row->len;
}

mr_float mr_sparse_row_get(mr_sparse_row *row, size_t col) {
    if (!row) {
        return 0.0f;
    }

    mr_index pos = MR_BINARY_SEARCH(row->cols, row->len, col);
    return pos != MR_INVALID_INDEX && row->cols[pos] == col ? row->values[pos] : 0.0f;
}

void mr_sparse_row_set(mr_sparse_row *row, size_t col, mr_float val) {
    if (!row) {
        return;
    }

    mr_index pos = MR_BINARY_SEARCH(row->cols, row->len, col);
    if (pos != MR_INVALID_INDEX && row->cols[pos] == col) {
        row->values[pos] = val;
    } else {
        row_insert(row, col, val, pos != MR_INVALID_INDEX ? pos : 0);
    }
}

void mr_sparse_row_clear(mr_sparse_row *row) {
    if (!row) {
        return;
    }

    row->len = 0;
}

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
    if (!builder) {
        return;
    }

    free(builder->rows);
    free(builder->cols);
    free(builder->values);

    free(builder);
}

void mr_sparse_matrix_builder_add_row(mr_sparse_matrix_builder *builder, mr_sparse_row *row) {
    if (!builder || !row || builder->nb_rows == builder->dim) {
        return;
    }

    size_t old_nb_cols = builder->nb_cols;
    builder->nb_cols += row->len;
    builder->rows[++builder->nb_rows] = builder->nb_cols;

    builder->values = xrealloc(builder->values, builder->nb_cols * sizeof(mr_float));
    memcpy(builder->values + old_nb_cols, row->values, row->len * sizeof(mr_float));

    builder->cols = xrealloc(builder->cols, builder->nb_cols * sizeof(size_t));
    memcpy(builder->cols + old_nb_cols, row->cols, row->len * sizeof(size_t));
}
