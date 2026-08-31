#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "maniray/utils/xmalloc.h"
#include "maniray/compute/math.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/compute/fvm/heat_dist.h"
#include "maniray/compute/fvm/scalar.h"
#include "maniray/compute/fvm/cell.h"

#define TIME_STEP 0.01

mr_fvm_heat_distance *mr_fvm_heat_distance_create() {
    mr_fvm_heat_distance *heat_distance = xmalloc(sizeof(mr_fvm_heat_distance));

    heat_distance->solver = mr_linear_system_solver_create();
    if (!heat_distance->solver) {
        free(heat_distance);

        return NULL;
    }

    heat_distance->forest = NULL;
    heat_distance->code_map = NULL;

    heat_distance->bc = NULL;
    heat_distance->init_cond_fn = NULL;
    heat_distance->sign_fn = NULL;

    heat_distance->heat_eq_mat = NULL;
    heat_distance->poisson_mat = NULL;

    heat_distance->init_cond_terms = NULL;

    return heat_distance;
}

void mr_fvm_heat_distance_destroy(mr_fvm_heat_distance *heat_distance) {
    if (!heat_distance) {
        return;
    }

    mr_vector_destroy(heat_distance->init_cond_terms);

    mr_sparse_matrix_destroy(heat_distance->poisson_mat);
    mr_sparse_matrix_destroy(heat_distance->heat_eq_mat);

    mr_boundary_condition_destroy(heat_distance->bc);

    mr_code_map_destroy(heat_distance->code_map);
    mr_ocforest_destroy(heat_distance->forest);

    mr_linear_system_solver_destroy(heat_distance->solver);

    free(heat_distance);
}

mr_ocforest *mr_fvm_heat_distance_ocforest_initialize(
    mr_fvm_heat_distance *heat_distance,
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots
) {
    assert(heat_distance);

    mr_code_map_destroy(heat_distance->code_map);
    heat_distance->code_map = NULL;

    mr_boundary_condition_destroy(heat_distance->bc);
    heat_distance->bc = NULL;

    mr_ocforest_destroy(heat_distance->forest);
    heat_distance->forest = mr_ocforest_create(
        manifold,
        roots,
        nb_roots,
        (size_t[]) { sizeof(mr_discretization_data), sizeof(mr_fvm_heat_distance_solution) },
        2
    );

    return heat_distance->forest;
}

mr_ocforest *mr_fvm_heat_distance_ocforest_update(mr_fvm_heat_distance *heat_distance) {
    assert(heat_distance);

    mr_code_map_destroy(heat_distance->code_map);
    heat_distance->code_map = NULL;

    return heat_distance->forest;
}

void mr_fvm_heat_distance_ocforest_finalize(mr_fvm_heat_distance *heat_distance) {
    assert(heat_distance);

    heat_distance->code_map = mr_code_map_create_from_ocforest(heat_distance->forest);
}

void mr_fvm_heat_distance_set_boundary_condition(mr_fvm_heat_distance *heat_distance, mr_boundary_condition *bc) {
    assert(heat_distance);

    mr_boundary_condition_destroy(heat_distance->bc);
    heat_distance->bc = bc;
}

void mr_fvm_heat_distance_set_initial_condition_fn(mr_fvm_heat_distance *heat_distance, mr_fvm_heat_distance_initial_cond_fn init_cond_fn) {
    assert(heat_distance);

    heat_distance->init_cond_fn = init_cond_fn;
}

void mr_fvm_heat_distance_set_sign_fn(mr_fvm_heat_distance *heat_distance, mr_fvm_heat_distance_sign_fn sign_fn) {
    assert(heat_distance);

    heat_distance->sign_fn = sign_fn;
}

typedef struct discr_matrix_data {
    mr_fvm_heat_distance *heat_distance;

    mr_sparse_row *temp_row;
    mr_sparse_matrix_builder *heat_eq_mat_builder;
    mr_sparse_matrix_builder *poisson_mat_builder;
} discr_matrix_data;

static void discr_matrix_data_create(mr_fvm_heat_distance *heat_distance, discr_matrix_data *data) {
    data->heat_distance = heat_distance;

    data->temp_row = mr_sparse_row_create();
    data->heat_eq_mat_builder = mr_sparse_matrix_builder_create(heat_distance->code_map->len);
    data->poisson_mat_builder = mr_sparse_matrix_builder_create(heat_distance->code_map->len);
}

static void discr_matrix_data_destroy(discr_matrix_data *data) {
    mr_sparse_matrix_builder_destroy(data->poisson_mat_builder);
    mr_sparse_matrix_builder_destroy(data->heat_eq_mat_builder);
    mr_sparse_row_destroy(data->temp_row);
}

static int write_matrix_coef(mr_ocforest *forest, mr_int cell_idx, mr_float64 coef, void *userdata) {
    discr_matrix_data *mat_data = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(mat_data->heat_distance->code_map, code);

    mr_float64 prev = mr_sparse_row_get(mat_data->temp_row, col);
    mr_sparse_row_set(mat_data->temp_row, col, prev + coef);

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
                    mat_data->heat_distance->bc,
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

    mr_sparse_matrix_builder_add_row(mat_data->poisson_mat_builder, mat_data->temp_row);

    mr_fvm_scalar_calc_transient_term(forest, cell_idx, TIME_STEP, store_cb, mr_fvm_scalar_store_coef_cb_null());
    mr_sparse_matrix_builder_add_row(mat_data->heat_eq_mat_builder, mat_data->temp_row);

    mr_sparse_row_clear(mat_data->temp_row);

    return MR_SUCCESS;
}

int mr_fvm_heat_distance_build_discretization_matrix(mr_fvm_heat_distance *heat_distance) {
    assert(heat_distance);

    discr_matrix_data mat_data;
    discr_matrix_data_create(heat_distance, &mat_data);

    for (mr_index octree_idx = 0; (size_t)octree_idx < heat_distance->forest->nb_roots; ++octree_idx) {
        int res = mr_octree_cells_apply(
            heat_distance->forest,
            octree_idx,
            mr_octree_apply_cb_create(fill_discr_matrix, &mat_data)
        );

        if (res != MR_SUCCESS) {
            discr_matrix_data_destroy(&mat_data);
            return res;
        }
    }

    heat_distance->heat_eq_mat = mr_sparse_matrix_build(mat_data.heat_eq_mat_builder);
    heat_distance->poisson_mat = mr_sparse_matrix_build(mat_data.poisson_mat_builder);

    mat_data.heat_eq_mat_builder = NULL;
    mat_data.poisson_mat_builder = NULL;
    discr_matrix_data_destroy(&mat_data);

    return MR_SUCCESS;
}

typedef struct init_cond_term_data {
    mr_fvm_heat_distance *heat_distance;

    mr_float64 *init_cond_term_arr;
} init_cond_term_data;

static int write_rhs_coef(mr_ocforest *forest, mr_int cell_idx, mr_float64 coef, void *userdata) {
    init_cond_term_data *init_cond_data = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(init_cond_data->heat_distance->code_map, code);

    init_cond_data->init_cond_term_arr[col] += coef;

    return MR_SUCCESS;
}

static int write_rhs_transient_coef(mr_ocforest *forest, mr_int cell_idx, mr_float64 coef, void *userdata) {
    init_cond_term_data *init_cond_data = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(init_cond_data->heat_distance->code_map, code);

    mr_float64 value = init_cond_data->heat_distance->init_cond_fn
        ? init_cond_data->heat_distance->init_cond_fn(init_cond_data->heat_distance, cell_idx)
        : 0.0f;
    init_cond_data->init_cond_term_arr[col] += coef * value;

    return MR_SUCCESS;
}

static int fill_source_term_array(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    init_cond_term_data *init_cond_data = userdata;

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type == MR_CELL_TYPE_EXTERIOR || discr_data->type == MR_CELL_TYPE_INTERPOLATION) {
        write_rhs_coef(forest, cell_idx, 0.0, init_cond_data);
    } else {
        mr_fvm_scalar_store_coef_cb store_cb = mr_fvm_scalar_store_coef_cb_create(write_rhs_coef, init_cond_data);
        for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
            if (!mr_is_boundary_cell(forest, cell_idx, dir)) {
                continue;
            }

            int res = mr_fvm_scalar_calc_boundary_flux(
                forest,
                init_cond_data->heat_distance->bc,
                cell_idx,
                dir,
                mr_fvm_scalar_store_coef_cb_null(),
                store_cb
            );

            if (res != MR_SUCCESS) {
                return res;
            }
        }

        mr_fvm_scalar_store_coef_cb store_transient_cb = mr_fvm_scalar_store_coef_cb_create(write_rhs_transient_coef, init_cond_data);
        mr_fvm_scalar_calc_transient_term(forest, cell_idx, TIME_STEP, mr_fvm_scalar_store_coef_cb_null(), store_transient_cb);
    }

    return MR_SUCCESS;
}

int mr_fvm_heat_distance_build_initial_condition_terms(mr_fvm_heat_distance *heat_distance) {
    assert(heat_distance);

    init_cond_term_data init_cond_data;
    init_cond_data.heat_distance = heat_distance;
    init_cond_data.init_cond_term_arr = xmalloc(heat_distance->code_map->len * sizeof(mr_float64));

    for (mr_index octree_idx = 0; (size_t)octree_idx < heat_distance->forest->nb_roots; ++octree_idx) {
        int res = mr_octree_cells_apply(
            heat_distance->forest,
            octree_idx,
            mr_octree_apply_cb_create(fill_source_term_array, &init_cond_data)
        );

        if (res != MR_SUCCESS) {
            free(init_cond_data.init_cond_term_arr);

            return res;
        }
    }

    heat_distance->init_cond_terms = mr_vector_create(init_cond_data.init_cond_term_arr, heat_distance->code_map->len);

    return MR_SUCCESS;
}

static int store_solution(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    mr_fvm_heat_distance *heat_distance = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(heat_distance->code_map, code);

    mr_fvm_heat_distance_solution *sol = mr_ocforest_get_cell_extra(forest, cell_idx, MR_HEAT_DIST_SOLUTION_EXTRA_FIELD);

    mr_float64 value = 0.0;
    int res = mr_linear_system_solver_get_solution(heat_distance->solver, col, &value);

    sol->dist = sqrt(4.0 * TIME_STEP * fabs(log(value)));

    return res;
}

int mr_fvm_heat_distance_solve(mr_fvm_heat_distance *heat_distance) {
    assert(heat_distance);

#if 0
    for (mr_int i = 0; i < heat_distance->heat_eq_mat->dim; ++i) {
        for (mr_int j = heat_distance->heat_eq_mat->rows[i]; j < heat_distance->heat_eq_mat->rows[i + 1]; ++j) {
            printf("%d:%f;", heat_distance->heat_eq_mat->cols[j], heat_distance->heat_eq_mat->values[j]);
        }
        printf("\n");
    }

    for (mr_int i = 0; i < heat_distance->init_cond_terms->len; ++i) {
        printf("%f\n", heat_distance->init_cond_terms->data[i]);
    }
#endif

    if (mr_linear_system_solver_set_matrix(heat_distance->solver, heat_distance->heat_eq_mat) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    if (mr_linear_system_solve(heat_distance->solver, heat_distance->init_cond_terms) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    for (mr_index octree_idx = 0; (size_t)octree_idx < heat_distance->forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(heat_distance->forest, octree_idx, mr_octree_apply_cb_create(store_solution, heat_distance));
    }

    return MR_SUCCESS;
}