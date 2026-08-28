#include <assert.h>
#include <stdio.h>

#include "maniray/utils/xmalloc.h"
#include "maniray/compute/math.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/compute/fvm/poisson.h"
#include "maniray/compute/fvm/scalar.h"
#include "maniray/compute/fvm/cell.h"

mr_fvm_poisson *mr_fvm_poisson_create() {
    mr_fvm_poisson *poisson = xmalloc(sizeof(mr_fvm_poisson));

    poisson->forest = NULL;
    poisson->code_map = NULL;

    poisson->source_fn = NULL;
    poisson->bc = NULL;

    poisson->discr_mat = NULL;
    poisson->source_terms = NULL;

    return poisson;
}

void mr_fvm_poisson_destroy(mr_fvm_poisson *poisson) {
    if (!poisson) {
        return;
    }

    mr_vector_destroy(poisson->source_terms);
    mr_sparse_matrix_destroy(poisson->discr_mat);

    mr_boundary_condition_destroy(poisson->bc);

    mr_code_map_destroy(poisson->code_map);
    mr_ocforest_destroy(poisson->forest);

    free(poisson);
}

mr_ocforest *mr_fvm_poisson_ocforest_initialize(
    mr_fvm_poisson *poisson,
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots
) {
    assert(poisson);

    mr_code_map_destroy(poisson->code_map);
    poisson->code_map = NULL;

    mr_boundary_condition_destroy(poisson->bc);
    poisson->bc = NULL;

    mr_ocforest_destroy(poisson->forest);
    poisson->forest = mr_ocforest_create(
        manifold,
        roots,
        nb_roots,
        (size_t[]) { sizeof(mr_discretization_data), sizeof(mr_fvm_poisson_solution) },
        2
    );

    return poisson->forest;
}

mr_ocforest *mr_fvm_poisson_ocforest_update(mr_fvm_poisson *poisson) {
    assert(poisson);

    mr_code_map_destroy(poisson->code_map);
    poisson->code_map = NULL;

    return poisson->forest;
}

void mr_fvm_poisson_ocforest_finalize(mr_fvm_poisson *poisson) {
    assert(poisson);

    poisson->code_map = mr_code_map_create_from_ocforest(poisson->forest);
}

void mr_fvm_poisson_set_source_term_fn(mr_fvm_poisson *poisson, mr_fvm_poisson_source_fn source_fn) {
    assert(poisson);

    poisson->source_fn = source_fn;
}

void mr_fvm_poisson_set_boundary_condition(mr_fvm_poisson *poisson, mr_boundary_condition *bc) {
    assert(poisson);

    mr_boundary_condition_destroy(poisson->bc);
    poisson->bc = bc;
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

static int write_matrix_coef(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    discr_matrix_data *mat_data = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    size_t col = mr_code_map_get_index(mat_data->poisson->code_map, code);

    mr_float prev = mr_sparse_row_get(mat_data->temp_row, col);
    // Subtract the coefficient so the matrix has the positive diagonal
    mr_sparse_row_set(mat_data->temp_row, col, prev - coef);

    return MR_SUCCESS;
}

static int fill_discr_matrix(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    discr_matrix_data *mat_data = userdata;

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    mr_fvm_scalar_store_coef_cb store_cb = mr_fvm_scalar_store_coef_cb_create(write_matrix_coef, mat_data);

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
                    mat_data->poisson->bc,
                    cell_idx,
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
    assert(poisson);

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

    poisson->discr_mat = mr_sparse_matrix_build(mat_data.matrix_builder);

    mat_data.matrix_builder = NULL;
    discr_matrix_data_destroy(&mat_data);

    return MR_SUCCESS;
}

typedef struct source_term_data {
    mr_fvm_poisson *poisson;

    mr_float64 *source_term_arr;
} source_term_data;

static int write_rhs_coef(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    source_term_data *src_data = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    size_t col = mr_code_map_get_index(src_data->poisson->code_map, code);

    // Subtract the coefficient to account for the matrix having a positive diagonal
    src_data->source_term_arr[col] -= coef;

    return MR_SUCCESS;
}

static int fill_source_term_array(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    source_term_data *src_data = userdata;

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type == MR_CELL_TYPE_EXTERIOR || discr_data->type == MR_CELL_TYPE_INTERPOLATION) {
        write_rhs_coef(forest, cell_idx, 0.0, src_data);
    } else {
        mr_fvm_scalar_store_coef_cb store_cb = mr_fvm_scalar_store_coef_cb_create(write_rhs_coef, src_data);
        for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
            if (!mr_is_boundary_cell(forest, cell_idx, dir)) {
                continue;
            }

            int res = mr_fvm_scalar_calc_boundary_flux(
                forest,
                src_data->poisson->bc,
                cell_idx,
                dir,
                mr_fvm_scalar_store_coef_cb_null(),
                store_cb
            );

            if (res != MR_SUCCESS) {
                return res;
            }
        }

        mr_float value = src_data->poisson->source_fn ? src_data->poisson->source_fn(src_data->poisson, cell_idx)
                                                      : 0.0f;
        mr_float volume = mr_cell_volume(forest, cell_idx);
        write_rhs_coef(forest, cell_idx, value * volume, src_data);
    }

    return MR_SUCCESS;
}

int mr_fvm_poisson_build_source_terms(mr_fvm_poisson *poisson) {
    assert(poisson);

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

    poisson->source_terms = mr_vector_create(st_data.source_term_arr, poisson->code_map->len);

    return MR_SUCCESS;
}

typedef struct store_solution_userdata {
    mr_fvm_poisson *poisson;
    LIS_VECTOR x;
} store_solution_userdata;

static int store_solution(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    store_solution_userdata *ud = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    size_t col = mr_code_map_get_index(ud->poisson->code_map, code);

    mr_fvm_poisson_solution *sol = mr_ocforest_get_cell_extra(forest, cell_idx, MR_POISSON_SOLUTION_EXTRA_FIELD);

    double value = 0.0;
    int res = lis_vector_get_value(ud->x, col, &value);
    sol->value = value;

    return res == LIS_SUCCESS ? MR_SUCCESS : MR_FAILURE;
}

int mr_fvm_poisson_solve(mr_fvm_poisson *poisson) {
    assert(poisson);

    LIS_VECTOR x;
    LIS_SOLVER solver;

    LIS_MATRIX A = poisson->discr_mat->inner;
    LIS_VECTOR b = poisson->source_terms->inner;

    lis_vector_duplicate(b, &x);
    lis_vector_copy(b, x);

    // Move to solver.c
    lis_solver_create(&solver);
    lis_solver_set_option("-initx_zeros 0 -i bicgstab -p ssor -tol 1.0e-8", solver);
    lis_solve(A, b, x, solver);

    store_solution_userdata ud = { poisson, x };
    for (mr_index octree_idx = 0; (size_t)octree_idx < poisson->forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(poisson->forest, octree_idx, mr_octree_apply_cb_create(store_solution, &ud));
    }

    lis_solver_destroy(solver);
    lis_vector_destroy(x);

    return MR_SUCCESS;
}