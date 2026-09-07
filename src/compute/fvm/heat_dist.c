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

// TODO: Make adaptive (equal to h^2)
#define TIME_STEP 0.0625

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

    if (cell_idx != 7679 || true) {
        mr_sparse_matrix_builder_add_row(mat_data->poisson_mat_builder, mat_data->temp_row);
    } else {
        mr_sparse_row *row = mr_sparse_row_create();

        mr_int code = mr_ocforest_get_code(forest, cell_idx);
        mr_int col = mr_code_map_get_index(mat_data->heat_distance->code_map, code);

        mr_sparse_row_set(row, col, 1.0);

        mr_sparse_matrix_builder_add_row(mat_data->poisson_mat_builder, row);
        mr_sparse_row_destroy(row);
    }

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

typedef struct store_gradient_component_userdata {
    mr_fvm_heat_distance *heat_distance;
    mr_float64 *grad_component;
} store_gradient_component_userdata;

static int store_gradient_component(mr_ocforest *forest, mr_int cell_idx, mr_float64 coef, void *userdata) {
    store_gradient_component_userdata *ud = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(ud->heat_distance->code_map, code);

    mr_float64 value = 0.0;
    if (mr_linear_system_solver_get_solution(ud->heat_distance->solver, col, &value) != MR_SUCCESS) {
        return MR_FAILURE;
    }
    *ud->grad_component += coef * value;

    return MR_SUCCESS;
}

typedef struct calculate_gradients_userdata {
    mr_fvm_heat_distance *heat_distance;
    mr_vector *poisson_rhs;
} calculate_gradients_userdata;

static int calculate_gradients(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    calculate_gradients_userdata *ud = userdata;

    mr_fvm_heat_distance *heat_distance = ud->heat_distance;
    mr_float64 grad[MR_NB_AXES] = { 0.0 };

    store_gradient_component_userdata store_ud[] = {
        { heat_distance, &grad[0] },
        { heat_distance, &grad[1] },
        { heat_distance, &grad[2] },
    };

    mr_fvm_scalar_store_coef_cb component_store[] = {
        mr_fvm_scalar_store_coef_cb_create(store_gradient_component, &store_ud[0]),
        mr_fvm_scalar_store_coef_cb_create(store_gradient_component, &store_ud[1]),
        mr_fvm_scalar_store_coef_cb_create(store_gradient_component, &store_ud[2]),
    };

    if (mr_fvm_scalar_calc_center_gradient(forest, cell_idx, component_store) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    /* mr_float64 grad_norm = mr_norm2(grad[0], grad[1], grad[2]);
    printf("%.16e\n", grad_norm); */

    /* if (grad_norm < 1.0e-12) {
        grad_norm = 1.0;
    } */

    mr_fvm_heat_distance_solution *sol = mr_ocforest_get_cell_extra(forest, cell_idx, MR_HEAT_DIST_SOLUTION_EXTRA_FIELD);
    // sol->dist = grad_norm;

    sol->grad[0] = grad[0];
    sol->grad[1] = grad[1];
    sol->grad[2] = grad[2];

#if 0
    sol->grad[0] = -grad[0] / grad_norm;
    sol->grad[1] = -grad[1] / grad_norm;
    sol->grad[2] = -grad[2] / grad_norm;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(heat_distance->code_map, code);

    // printf("%.16f, %.16f, %.16f\n", ud->poisson_rhs->data[col] * 100000.0, grad_norm * 100000.0, -ud->poisson_rhs->data[col] / grad_norm);

    if (cell_idx != 7679) {
        ud->poisson_rhs->data[col] = -ud->poisson_rhs->data[col] / grad_norm;
    } else {
        ud->poisson_rhs->data[col] = 0.0;

        printf("CELL %d -- %d\n", cell_idx, col);
    }

    printf("%f\n", ud->poisson_rhs->data[col]);


    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    if (mr_is_boundary_cell(forest, cell_idx, MR_DIRECTION_PL_X) && MR_ABS(cell->y) < 0.5f && MR_ABS(cell->z) < 0.5f) {
        mr_float64 m = 1e18;
        printf("Cell %d: (%f, %f, %f)\n", cell_idx, grad[0] * m, grad[1] * m, grad[2] * m);
    }
#endif

    return MR_SUCCESS;
}

typedef struct calc_flux_userdata {
    mr_fvm_heat_distance *heat_distance;
    mr_float64 value;
} calc_flux_userdata;

static int calc_flux(mr_ocforest *forest, mr_int cell_idx, mr_float64 coef, void *userdata) {
    calc_flux_userdata *ud = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(ud->heat_distance->code_map, code);

    mr_float64 value = 0.0;
    if (mr_linear_system_solver_get_solution(ud->heat_distance->solver, col, &value) != MR_SUCCESS) {
        return MR_FAILURE;
    }
    ud->value += value * coef;

    return MR_SUCCESS;
}

typedef struct calc_face_gradient_norm_userdata {
    mr_fvm_heat_distance *heat_distance;
    mr_float64 value;
    mr_float64 grad[MR_NB_AXES];
} calc_face_gradient_norm_userdata;

static int calc_face_gradient_norm(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    calc_face_gradient_norm_userdata *ud = userdata;

    mr_fvm_heat_distance_solution *sol = mr_ocforest_get_cell_extra(forest, cell_idx, MR_HEAT_DIST_SOLUTION_EXTRA_FIELD);
    mr_float64 grad_norm = mr_norm2(sol->grad[0], sol->grad[1], sol->grad[2]);

    ud->value += grad_norm * coef;
    // printf("%d: %f\n", cell_idx, coef);

    ud->grad[0] += sol->grad[0] * coef;
    ud->grad[1] += sol->grad[1] * coef;
    ud->grad[2] += sol->grad[2] * coef;

    return MR_SUCCESS;
}

typedef struct fill_poisson_rhs_userdata {
    mr_fvm_heat_distance *heat_distance;
    mr_vector *poisson_rhs;
} fill_poisson_rhs_userdata;

static int fill_poisson_rhs(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    fill_poisson_rhs_userdata *ud = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(ud->heat_distance->code_map, code);

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

    if (discr_data->type == MR_CELL_TYPE_EXTERIOR || discr_data->type == MR_CELL_TYPE_INTERPOLATION /* || cell_idx == 7679 */) {
        ud->poisson_rhs->data[col] = 0.0;
    } else {
        for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
            calc_flux_userdata flux_store_ud = { ud->heat_distance, 0.0 };
            mr_fvm_scalar_store_coef_cb flux_store_cb = mr_fvm_scalar_store_coef_cb_create(calc_flux, &flux_store_ud);

            int res = MR_SUCCESS;
            if (mr_is_boundary_cell(forest, cell_idx, dir)) {
                res = mr_fvm_scalar_calc_boundary_flux(
                    forest,
                    ud->heat_distance->bc,
                    cell_idx,
                    dir,
                    flux_store_cb,
                    mr_fvm_scalar_store_coef_cb_null() // TODO: Add store for it
                );
            } else {
                res = mr_fvm_scalar_calc_internal_flux(forest, cell_idx, dir, flux_store_cb);
            }

            if (res != MR_SUCCESS) {
                return res;
            }

            calc_face_gradient_norm_userdata face_grad_store_ud = { ud->heat_distance, 0.0, { 0.0 } };
            mr_fvm_interpolation_cb face_grad_store_cb = mr_fvm_interpolation_cb_create(calc_face_gradient_norm, &face_grad_store_ud);

            if (mr_fvm_interpolate_face_value(forest, cell_idx, dir, face_grad_store_cb) != MR_SUCCESS) {
                return MR_FAILURE;
            }

            // mr_float64 norm = mr_norm2(face_grad_store_ud.grad[0], face_grad_store_ud.grad[1], face_grad_store_ud.grad[2]);
            ud->poisson_rhs->data[col] += -flux_store_ud.value / face_grad_store_ud.value;
        }
    }

    // printf("%.16e\n", ud->poisson_rhs->data[col]);

    if (0) {
        mr_float64 value = 0.0;
        mr_linear_system_solver_get_solution(ud->heat_distance->solver, col, &value);
        printf("%.16e\n", value);
    }

    return MR_SUCCESS;
}

static int store_solution(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    mr_fvm_heat_distance *heat_distance = userdata;

    mr_int code = mr_ocforest_get_code(forest, cell_idx);
    mr_int col = mr_code_map_get_index(heat_distance->code_map, code);

    // mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_fvm_heat_distance_solution *sol = mr_ocforest_get_cell_extra(forest, cell_idx, MR_HEAT_DIST_SOLUTION_EXTRA_FIELD);

    mr_float64 value = 0.0;
    int res = mr_linear_system_solver_get_solution(heat_distance->solver, col, &value);
    sol->dist = value;

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

    mr_int poisson_rhs_len = heat_distance->init_cond_terms->len;
    mr_float64 *poisson_rhs_data = xcalloc(poisson_rhs_len, sizeof(mr_float64));
    mr_vector *poisson_rhs = mr_vector_create(poisson_rhs_data, poisson_rhs_len);

    calculate_gradients_userdata calc_grad_ud = { heat_distance, poisson_rhs };
    for (mr_index octree_idx = 0; (size_t)octree_idx < heat_distance->forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(heat_distance->forest, octree_idx, mr_octree_apply_cb_create(calculate_gradients, &calc_grad_ud));
    }

    fill_poisson_rhs_userdata fill_rhs_ud = { heat_distance, poisson_rhs };
    for (mr_index octree_idx = 0; (size_t)octree_idx < heat_distance->forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(heat_distance->forest, octree_idx, mr_octree_apply_cb_create(fill_poisson_rhs, &fill_rhs_ud));
    }

    if (mr_linear_system_solver_set_matrix(heat_distance->solver, heat_distance->poisson_mat) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    if (mr_linear_system_solve(heat_distance->solver, poisson_rhs) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    for (mr_index octree_idx = 0; (size_t)octree_idx < heat_distance->forest->nb_roots; ++octree_idx) {
        mr_octree_cells_apply(heat_distance->forest, octree_idx, mr_octree_apply_cb_create(store_solution, heat_distance));
    }

    mr_vector_destroy(poisson_rhs);

    return MR_SUCCESS;
}