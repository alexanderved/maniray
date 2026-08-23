#include <stdio.h>
#include <tgmath.h>

#include "maniray/utils/xmalloc.h"
#include "maniray/compute/fvm/boundary.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/compute/fvm/poisson.h"
#include "maniray/compute/fvm/scalar.h"
#include "maniray/compute/matrix.h"

mr_fvm_poisson *mr_fvm_poisson_create(
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots,
    mr_fvm_poisson_source_fn source_fn
) {
    mr_fvm_poisson *poisson = xmalloc(sizeof(mr_fvm_poisson));

    poisson->source_fn = source_fn;

    poisson->forest = mr_ocforest_create(
        manifold,
        roots,
        nb_roots,
        (size_t[]) { sizeof(mr_discretization_data), sizeof(mr_fvm_poisson_solution) },
        2
    );
    poisson->code_map = NULL;

    poisson->discr_mat = NULL;
    poisson->source_terms = NULL;
    poisson->res = NULL;

    return poisson;
}

void mr_fvm_poisson_destroy(mr_fvm_poisson *poisson) {
    if (!poisson) {
        return;
    }

    free(poisson->res);
    mr_dense_matrix_destroy(poisson->source_terms);
    mr_sparse_matrix_destroy(poisson->discr_mat);

    mr_code_map_destroy(poisson->code_map);
    mr_ocforest_destroy(poisson->forest);

    free(poisson);
}

mr_ocforest *mr_fvm_poisson_get_ocforest(mr_fvm_poisson *poisson) {
    return poisson ? poisson->forest : NULL;
}

typedef struct discr_matrix_data {
    mr_fvm_poisson *poisson;

    mr_sparse_row *temp_row;
    mr_sparse_matrix_builder *matrix_builder;
} discr_matrix_data;

static void discr_matrix_data_create(mr_fvm_poisson *poisson, discr_matrix_data *data) {
    data->poisson = poisson;

    data->temp_row = mr_sparse_row_create();
    data->matrix_builder = mr_sparse_matrix_builder_create(poisson->code_map->len);
}

static void discr_matrix_data_destroy(discr_matrix_data *data) {
    mr_sparse_matrix_builder_destroy(data->matrix_builder);
    mr_sparse_row_destroy(data->temp_row);
}

static int write_coef(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    discr_matrix_data *mat_data = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    size_t col = mr_code_map_get_index(mat_data->poisson->code_map, code);

    // printf("%d: %d\n", cell_idx, code);

    mr_float prev = mr_sparse_row_get(mat_data->temp_row, col);
    mr_sparse_row_set(mat_data->temp_row, col, prev + coef);

    return MR_SUCCESS;
}

static int fill_discr_matrix(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    discr_matrix_data *mat_data = userdata;

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    mr_fvm_scalar_store_coef_cb store_cb = mr_fvm_scalar_store_coef_cb_create(write_coef, mat_data);

    int res = MR_SUCCESS;
    if (discr_data->type == MR_CELL_TYPE_EXTERIOR) {
        res = mr_fvm_scalar_mark_inactive_cell(forest, cell_idx, store_cb);
    } else if (discr_data->type == MR_CELL_TYPE_INTERPOLATION) {
        res = mr_fvm_scalar_interpolate(forest, cell_idx, store_cb);
    } else {
        for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
            if (mr_is_boundary_cell(forest, cell_idx, dir)) {
                res = mr_fvm_scalar_calc_boundary_flux(
                    forest,
                    cell_idx,
                    NULL,
                    dir,
                    store_cb,
                    mr_fvm_scalar_store_coef_cb_null()
                );
            } else {
                res = mr_fvm_scalar_calc_internal_flux(forest, cell_idx, dir, store_cb);
            }

            if (res != MR_SUCCESS) {
                return res;
            }
        }
    }

    if (res != MR_SUCCESS) {
        return res;
    }

    mr_sparse_matrix_builder_add_row(mat_data->matrix_builder, mat_data->temp_row);
    mr_sparse_row_clear(mat_data->temp_row);

    return MR_SUCCESS;
}

int mr_fvm_poisson_build_discretization_matrix(mr_fvm_poisson *poisson) {
    if (!poisson) {
        return MR_FAILURE;
    }

    // TODO: create functions to start and end building octree and move this code there
    if (!poisson->code_map) {
        poisson->code_map = mr_code_map_create_from_ocforest(poisson->forest);
    }

    discr_matrix_data mat_data;
    discr_matrix_data_create(poisson, &mat_data);

    for (mr_index octree_idx = 0; (size_t)octree_idx < poisson->forest->nb_roots; ++octree_idx) {
        int res = mr_octree_cells_apply(
            poisson->forest,
            octree_idx,
            mr_octree_apply_cb_create(fill_discr_matrix, &mat_data)
        );

        if (res != MR_SUCCESS) {
            discr_matrix_data_destroy(&mat_data);
            return res;
        }
    }

#if 0
    for (mr_int i = 0; i < mat_data.matrix_builder->dim; ++i) {
        printf("Row %d: ", i);
        for (mr_int j = mat_data.matrix_builder->rows[i]; j < mat_data.matrix_builder->rows[i + 1]; ++j) {
            printf("%d: %f\t", mat_data.matrix_builder->cols[j], mat_data.matrix_builder->values[j]);
        }
        printf("\n");
    }
#endif

    poisson->discr_mat = mr_sparse_matrix_build(mat_data.matrix_builder);

    mat_data.matrix_builder = NULL;
    discr_matrix_data_destroy(&mat_data);

    return MR_SUCCESS;
}

typedef struct source_term_data {
    mr_fvm_poisson *poisson;

    mr_float64 *source_term_arr;
} source_term_data;

// TODO: Add boundary handling
static int fill_source_term_array(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    source_term_data *src_data = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    size_t col = mr_code_map_get_index(src_data->poisson->code_map, code);

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type == MR_CELL_TYPE_EXTERIOR || discr_data->type == MR_CELL_TYPE_INTERPOLATION) {
        src_data->source_term_arr[col] = 0.0f;
    } else {
        mr_float value = src_data->poisson->source_fn ? src_data->poisson->source_fn(src_data->poisson, cell_idx) : 0.0f;
        src_data->source_term_arr[col] = value;
    }

    return MR_SUCCESS;
}

int mr_fvm_poisson_build_source_terms(mr_fvm_poisson *poisson) {
    source_term_data st_data;
    st_data.poisson = poisson;
    st_data.source_term_arr = xmalloc(poisson->code_map->len * sizeof(mr_float64));

    for (mr_index octree_idx = 0; (size_t)octree_idx < poisson->forest->nb_roots; ++octree_idx) {
        int res = mr_octree_cells_apply(
            poisson->forest,
            octree_idx,
            mr_octree_apply_cb_create(fill_source_term_array, &st_data)
        );

        if (res != MR_SUCCESS) {
            free(st_data.source_term_arr);

            return res;
        }
    }

    poisson->source_terms = mr_dense_matrix_create(st_data.source_term_arr, poisson->code_map->len);

    return MR_SUCCESS;
}

int mr_fvm_poisson_solve(mr_fvm_poisson *poisson) {
    // TODO: Handle errors better
    Vec x;
    KSP ksp;

    Mat A = poisson->discr_mat->mat;
    Vec b = poisson->source_terms->vec;

    if (VecDuplicate(b, &x)) abort();
    if (VecSet(x, 0.0)) abort();

    if (KSPCreate(PETSC_COMM_SELF, &ksp)) abort();

    // Required to make the matrix non-singular
    MatNullSpace nsp;
    MatNullSpaceCreate(PETSC_COMM_SELF, PETSC_TRUE, 0, NULL, &nsp);
    MatSetNullSpace(A, nsp);
    MatSetTransposeNullSpace(A, nsp);

    if (KSPSetOperators(ksp, A, A)) abort();

    if (KSPSetType(ksp, KSPGMRES)) abort();

    /* PC pc;
    if (KSPGetPC(ksp, &pc)) abort();
    if (PCSetType(pc, PCGAMG)) abort(); */

    if (KSPSetFromOptions(ksp)) abort();
    if (KSPSolve(ksp, b, x)) abort();

    const mr_float64 *res_arr;
    if (VecGetArrayRead(x, &res_arr)) abort();

    poisson->res = xmalloc(poisson->source_terms->len * sizeof(mr_float64));
    memcpy(poisson->res, res_arr, poisson->source_terms->len * sizeof(mr_float64));

    if (VecRestoreArrayRead(x, &res_arr)) abort();

    // if (VecView(x, PETSC_VIEWER_STDOUT_SELF)) abort();


    MatNullSpaceDestroy(&nsp);
    if (KSPDestroy(&ksp)) abort();
    if (VecDestroy(&x)) abort();


#if 0
    for (mr_int i = 0; i < poisson->discr_mat->dim; ++i) {
        printf("%d: %f\n", i, poisson->res[i]);
    }
#endif

    return MR_SUCCESS;
}