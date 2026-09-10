#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "maniray/compute/fvm/scalar.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/geometry.h"

typedef struct store_with_mul_userdata {
    mr_fvm_scalar_store_coef_cb store;
    mr_float64 mul;
} store_with_mul_userdata;

static int store_with_mul(mr_ocforest *forest, mr_int cell_idx, mr_float64 coef, void *userdata) {
    store_with_mul_userdata *ud = userdata;

    return ud->store.fn(forest, cell_idx, coef * ud->mul, ud->store.userdata);
}

static int calc_forward_diff_derivative(
    mr_ocforest *forest,
    mr_int part_stencil[2],
    mr_axis axis,
    mr_fvm_scalar_store_coef_cb store
) {
    mr_int full_stencil[3] = { part_stencil[0], part_stencil[1], MR_INVALID_INDEX };

    mr_octree_cell_neighbor neighbor = mr_octree_find_face_neighbor_cells(
        forest,
        full_stencil[1],
        mr_direction_create(axis, MR_SIGN_PLUS)
    );
    assert(neighbor.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);
    full_stencil[2] = neighbor.neighbor_idx;

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, full_stencil[0]);
    mr_float64 divisor = 1.0 / cell->dim;

    store.fn(forest, full_stencil[0], -1.5 * divisor, store.userdata);
    store.fn(forest, full_stencil[1], 2.0 * divisor, store.userdata);
    store.fn(forest, full_stencil[2], -0.5 * divisor, store.userdata);

    return MR_SUCCESS;
}

static int calc_backward_diff_derivative(
    mr_ocforest *forest,
    mr_int part_stencil[2],
    mr_axis axis,
    mr_fvm_scalar_store_coef_cb store
) {
    mr_int full_stencil[3] = { MR_INVALID_INDEX, part_stencil[0], part_stencil[1] };

    mr_octree_cell_neighbor neighbor = mr_octree_find_face_neighbor_cells(
        forest,
        full_stencil[1],
        mr_direction_create(axis, MR_SIGN_MINUS)
    );
    assert(neighbor.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);
    full_stencil[0] = neighbor.neighbor_idx;

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, full_stencil[2]);
    mr_float64 divisor = 1.0 / cell->dim;

    store.fn(forest, full_stencil[0], 0.5 * divisor, store.userdata);
    store.fn(forest, full_stencil[1], -2.0 * divisor, store.userdata);
    store.fn(forest, full_stencil[2], 1.5 * divisor, store.userdata);

    return MR_SUCCESS;
}

static int calc_center_diff_derivative(
    mr_ocforest *forest,
    mr_int stencil[2],
    mr_fvm_scalar_store_coef_cb store
) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, stencil[0]);
    mr_float64 divisor = 1.0 / (2.0 * cell->dim);

    store.fn(forest, stencil[0], -divisor, store.userdata);
    store.fn(forest, stencil[1], divisor, store.userdata);

    return MR_SUCCESS;
}

int mr_fvm_scalar_calc_center_derivative(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_axis axis,
    mr_fvm_scalar_store_coef_cb store
) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);
    assert(!mr_fvm_scalar_store_coef_cb_is_null(store));

    mr_octree_cell_neighbor backward_neighbor = mr_octree_find_face_neighbor_cells(
        forest,
        cell_idx,
        mr_direction_create(axis, MR_SIGN_MINUS)
    );
    mr_octree_cell_neighbor forward_neighbor = mr_octree_find_face_neighbor_cells(
        forest,
        cell_idx,
        mr_direction_create(axis, MR_SIGN_PLUS)
    );

    if (backward_neighbor.type != MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE
        || !mr_ocforest_is_node_active(forest, backward_neighbor.node_idx))
    {
        assert(forward_neighbor.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);
        return calc_forward_diff_derivative(
            forest,
            (mr_int[]) { cell_idx, forward_neighbor.neighbor_idx },
            axis,
            store
        );
    }

    if (forward_neighbor.type != MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE
        || !mr_ocforest_is_node_active(forest, forward_neighbor.node_idx))
    {
        assert(backward_neighbor.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);
        return calc_backward_diff_derivative(
            forest,
            (mr_int[]) { backward_neighbor.neighbor_idx, cell_idx },
            axis,
            store
        );
    }

    assert(backward_neighbor.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);
    assert(forward_neighbor.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);

    return calc_center_diff_derivative(
        forest,
        (mr_int[]) { backward_neighbor.neighbor_idx, forward_neighbor.neighbor_idx },
        store
    );
}

typedef struct interpolate_derivative_to_face_userdata {
    mr_fvm_scalar_store_coef_cb store;
    mr_axis axis;
    mr_float64 mul;
} interpolate_derivative_to_face_userdata;

static int interpolate_derivative_to_face(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    interpolate_derivative_to_face_userdata *ud = userdata;

    return mr_fvm_scalar_calc_center_derivative(
        forest,
        cell_idx,
        ud->axis,
        mr_fvm_scalar_store_coef_cb_create(store_with_mul, &(store_with_mul_userdata) { ud->store, coef * ud->mul })
    );
}

static int interpolate_transverse_derivative(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_direction face_dir,
    mr_float mul,
    mr_fvm_scalar_store_coef_cb store
) {
    interpolate_derivative_to_face_userdata ud = { store, mr_direction_get_axis(face_dir), mul };
    mr_fvm_interpolation_cb interp = mr_fvm_interpolation_cb_create(interpolate_derivative_to_face, &ud);
    
    return mr_fvm_interpolate_face_value(forest, cell_idx, face_dir, interp);
}

static int calc_equal_size_transverse_derivative(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_int neighbor_cell_idx,
    mr_direction dir,
    mr_float mul,
    mr_fvm_scalar_store_coef_cb store
) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_float coef = mul / cell->dim * mr_direction_get_sign_mul(dir);

    if (store.fn(forest, cell_idx, -coef, store.userdata) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    return store.fn(forest, neighbor_cell_idx, coef, store.userdata);
}

static int store_transverse_der_ghost_cell_coefs(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    store_with_mul_userdata *gc_ud = userdata;
    return gc_ud->store.fn(forest, cell_idx, coef * gc_ud->mul, gc_ud->store.userdata);
}

static int calc_fine_coarse_transverse_derivative(
    mr_ocforest *forest,
    mr_int fine_cell_idx,
    mr_int coarse_cell_idx,
    mr_direction dir,
    mr_float mul,
    mr_fvm_scalar_store_coef_cb store
) {
    mr_octree_cell *fine_cell = mr_ocforest_get_cell(forest, fine_cell_idx);
    mr_float coef = mul / fine_cell->dim * mr_direction_get_sign_mul(dir);

    if (store.fn(forest, fine_cell_idx, -coef, store.userdata) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    store_with_mul_userdata ud = { store, coef };
    return mr_fvm_calculate_ghost_cell(
        forest,
        coarse_cell_idx,
        fine_cell_idx,
        mr_fvm_interpolation_cb_create(store_transverse_der_ghost_cell_coefs, &ud)
    );
}

static int calc_transverse_derivative(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_direction face_dir,
    mr_fvm_scalar_store_coef_cb store
) {
    mr_octree_cell_neighbor cell_neighbor = mr_octree_find_face_neighbor_cells(forest, cell_idx, face_dir);
    if (cell_neighbor.node_idx == MR_INVALID_INDEX || !mr_ocforest_is_node_active(forest, cell_neighbor.node_idx)) {
        return interpolate_transverse_derivative(forest, cell_idx, face_dir, 1.0f, store);
    }

    switch (cell_neighbor.type) {
        case MR_OCTREE_CELL_NEIGHBOR_NONE:
            return MR_FAILURE;

        case MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE:
            return calc_equal_size_transverse_derivative(
                forest,
                cell_idx,
                cell_neighbor.neighbor_idx,
                face_dir,
                1.0f,
                store
            );

        case MR_OCTREE_CELL_NEIGHBOR_COARSER:
            return calc_fine_coarse_transverse_derivative(
                forest,
                cell_idx,
                cell_neighbor.neighbor_idx,
                face_dir,
                1.0f,
                store
            );

        case MR_OCTREE_CELL_NEIGHBOR_FINER: ;
            mr_direction reflected_dir = mr_direction_reflect(face_dir);

            mr_float fine_face_areas[4] = { 0.0f };
            mr_float coarse_face_area = 0.0f;

            for (size_t i = 0; i < 4; ++i) {
                fine_face_areas[i] = mr_cell_face_area(forest, cell_neighbor.neighbor_indices[i], reflected_dir);
                coarse_face_area += fine_face_areas[i];
            }

            for (size_t i = 0; i < 4; ++i) {
                if (calc_fine_coarse_transverse_derivative(
                    forest,
                    cell_neighbor.neighbor_indices[i],
                    cell_idx,
                    reflected_dir,
                    fine_face_areas[i] / coarse_face_area,
                    store
                ) != MR_SUCCESS) {
                    return MR_FAILURE;
                }
            }

            return MR_SUCCESS;
    }

    return MR_FAILURE;
}

static int calc_tangent_derivative(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_direction face_dir,
    mr_axis axis,
    mr_fvm_scalar_store_coef_cb store
) {
    interpolate_derivative_to_face_userdata ud = { store, axis, 1.0 };
    mr_fvm_interpolation_cb interp = mr_fvm_interpolation_cb_create(interpolate_derivative_to_face, &ud);
    
    return mr_fvm_interpolate_face_value(forest, cell_idx, face_dir, interp);
}

int mr_fvm_scalar_calc_face_derivative(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_direction face_dir,
    mr_axis axis,
    mr_fvm_scalar_store_coef_cb store
) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);
    assert(!mr_fvm_scalar_store_coef_cb_is_null(store));

    if (mr_direction_get_axis(face_dir) == axis) {
        return calc_transverse_derivative(forest, cell_idx, face_dir, store);
    }

    return calc_tangent_derivative(forest, cell_idx, face_dir, axis, store);
}

static void precompute_inv_metric(mr_manifold *manifold, mr_uint chart_idx, mr_float p[MR_NB_AXES], mr_float inv_metric[6]) {
    for (size_t i = 0; i < MR_NB_AXES; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            inv_metric[i * (i + 1) / 2 + j] = mr_manifold_inv_metric(manifold, chart_idx, p, i, j);
        }
    }
}

static mr_float get_from_precomputed_inv_metric(mr_float inv_metric[6], size_t i, size_t j) {
    if (i < j) {
        size_t tmp = i;
        i = j;
        j = tmp;
    }

    return inv_metric[i * (i + 1) / 2 + j];
}

typedef struct store_derivative_in_grad_userdata {
    mr_fvm_scalar_store_coef_cb *component_store;
    mr_axis axis;
    mr_float inv_metric[6];
} store_derivative_in_grad_userdata;

static int store_derivative_in_grad(mr_ocforest *forest, mr_int cell_idx, mr_float64 coef, void *userdata) {
    store_derivative_in_grad_userdata *ud = userdata;
    mr_fvm_scalar_store_coef_cb *component_store = ud->component_store;

    for (mr_int i = 0; i < MR_NB_AXES; ++i) {
        mr_float inv_metric_comp = get_from_precomputed_inv_metric(ud->inv_metric, i, ud->axis);
        if (component_store[i].fn(forest, cell_idx, inv_metric_comp * coef, component_store[i].userdata) != MR_SUCCESS) {
            return MR_FAILURE;
        }
    }

    return MR_SUCCESS;
}

int mr_fvm_scalar_calc_center_gradient(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_fvm_scalar_store_coef_cb component_store[MR_NB_AXES]
) {
    assert(forest);
    assert(component_store);
    assert(cell_idx != MR_INVALID_INDEX);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_float center[] = { cell->x, cell->y, cell->z };
    
    store_derivative_in_grad_userdata ud = {
        .component_store = component_store,
    };
    precompute_inv_metric(forest->manifold, cell->chart_idx, center, ud.inv_metric);

    mr_fvm_scalar_store_coef_cb der_store = mr_fvm_scalar_store_coef_cb_create(store_derivative_in_grad, &ud);
    for (mr_int j = 0; j < MR_NB_AXES; ++j) {
        ud.axis = j;
        if (mr_fvm_scalar_calc_center_derivative(forest, cell_idx, j, der_store) != MR_SUCCESS) {
            return MR_FAILURE;
        }
    }

    return MR_SUCCESS;
}

int mr_fvm_scalar_calc_face_gradient(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_direction face_dir,
    mr_fvm_scalar_store_coef_cb component_store[MR_NB_AXES]
) {
    assert(forest);
    assert(component_store);
    assert(cell_idx != MR_INVALID_INDEX);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

    mr_float face_center[MR_NB_AXES] = { 0.0f };
    mr_cell_face_center(forest, cell_idx, face_dir, face_center);
    
    store_derivative_in_grad_userdata ud = {
        .component_store = component_store,
    };
    precompute_inv_metric(forest->manifold, cell->chart_idx, face_center, ud.inv_metric);

    mr_fvm_scalar_store_coef_cb der_store = mr_fvm_scalar_store_coef_cb_create(store_derivative_in_grad, &ud);
    for (mr_int j = 0; j < MR_NB_AXES; ++j) {
        ud.axis = j;
        if (mr_fvm_scalar_calc_face_derivative(forest, cell_idx, face_dir, j, der_store) != MR_SUCCESS) {
            return MR_FAILURE;
        }
    }

    return MR_SUCCESS;
}

int mr_fvm_scalar_mark_inactive_cell(mr_ocforest *forest, mr_int cell_idx, mr_fvm_scalar_store_coef_cb store) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);

    if (mr_fvm_scalar_store_coef_cb_is_null(store)) {
        return MR_SUCCESS;
    }

    return store.fn(forest, cell_idx, 1.0f, store.userdata);
}

static int get_interp_point(mr_ocforest *forest, mr_int cell_idx, mr_index *octree_idx, mr_float p[MR_NB_AXES]) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);

    *octree_idx = discr_data->donor_root_idx;
    mr_int donor_root_node_idx = mr_ocforest_get_root_node(forest, discr_data->donor_root_idx);
    mr_uint donor_chart = mr_ocforest_get_node(forest, donor_root_node_idx)->chart_idx;

    mr_float center[MR_NB_AXES] = { cell->x, cell->y, cell->z };
    return mr_manifold_transition(forest->manifold, cell->chart_idx, donor_chart, p, center);
}

static int store_interpolation_coefs(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    mr_fvm_scalar_store_coef_cb *store = userdata;

    return store->fn(forest, cell_idx, -coef, store->userdata);
}

int mr_fvm_scalar_interpolate(mr_ocforest *forest, mr_int cell_idx, mr_fvm_scalar_store_coef_cb store) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);

    if (mr_fvm_scalar_store_coef_cb_is_null(store)) {
        return MR_SUCCESS;
    }

    if (store.fn(forest, cell_idx, 1.0f, store.userdata) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    mr_index octree_idx = MR_INVALID_INDEX;
    mr_float p[MR_NB_AXES] = { 0.0f };
    if (get_interp_point(forest, cell_idx, &octree_idx, p) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    return mr_fvm_perform_interpolation(forest, octree_idx, p, mr_fvm_interpolation_cb_create(store_interpolation_coefs, &store));
}

int mr_fvm_scalar_calc_boundary_flux(
    mr_ocforest *forest,
    mr_boundary_condition *bc,
    mr_int cell_idx,
    mr_direction dir,
    mr_fvm_scalar_store_coef_cb store_implicit,
    mr_fvm_scalar_store_coef_cb store_rhs
) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);
    assert(bc);

    switch (mr_boundary_condition_get_type(bc, cell_idx, dir)) {
        // Direchlet BC is implemented using second-order one-sided finite differences
        case MR_BC_DIRICHLET: ;
            mr_axis axis = mr_direction_get_axis(dir);
            mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

            mr_float face_center[MR_NB_AXES] = { 0.0f };
            mr_cell_face_center(forest, cell_idx, dir, face_center);

            mr_float sqrt_inv_coef = sqrtf(mr_manifold_inv_metric(forest->manifold, cell->chart_idx, face_center, axis, axis));
            mr_float area = mr_cell_face_area(forest, cell_idx, dir);

            if (!mr_fvm_scalar_store_coef_cb_is_null(store_implicit)) {
                mr_octree_cell_neighbor cell_neighbor = mr_octree_find_face_neighbor_cells(forest, cell_idx, mr_direction_reflect(dir));
                assert(cell_neighbor.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);

                // TODO: Handle non-orthogonal curvilinear coordinates
                mr_float cell_coef = -sqrt_inv_coef * 3.0f / cell->dim * area;
                mr_float neighbor_coef = sqrt_inv_coef * 1.0f / (3.0f * cell->dim) * area;

                store_implicit.fn(forest, cell_idx, -cell_coef, store_implicit.userdata);
                store_implicit.fn(forest, cell_neighbor.neighbor_idx, -neighbor_coef, store_implicit.userdata);
            }

            if (!mr_fvm_scalar_store_coef_cb_is_null(store_rhs)) {
                mr_float value = 0.0f;
                mr_boundary_condition_get_value(bc, cell_idx, dir, &value);

                mr_float coef = sqrt_inv_coef * 8.0f * value / (3.0f * cell->dim) * area;
                store_rhs.fn(forest, cell_idx, coef, store_rhs.userdata);
            }

            break;

        case MR_BC_NEUMANN:
            if (!mr_fvm_scalar_store_coef_cb_is_null(store_rhs)) {
                mr_float area = mr_cell_face_area(forest, cell_idx, dir);
                mr_float value = 0.0f;
                mr_boundary_condition_get_value(bc, cell_idx, dir, &value);

                mr_float coef = value * area;
                store_rhs.fn(forest, cell_idx, coef, store_rhs.userdata);
            }

            break;

        default:
            assert(false);
    }

    return MR_SUCCESS;
}

static mr_float calc_transverse_der_mul(mr_ocforest *forest, mr_int cell_idx, mr_direction dir) {
    mr_axis axis = mr_direction_get_axis(dir);
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

    mr_float middle[MR_NB_AXES] = { 0.0f };
    mr_cell_face_center(forest, cell_idx, dir, middle);

    mr_float sqrt_inv_coef = sqrt(mr_manifold_inv_metric(forest->manifold, cell->chart_idx, middle, axis, axis));
    mr_float area = mr_cell_face_area(forest, cell_idx, dir);

    return -sqrt_inv_coef * area * mr_direction_get_sign_mul(dir);
}

int mr_fvm_scalar_calc_internal_flux(mr_ocforest *forest, mr_int cell_idx, mr_direction dir, mr_fvm_scalar_store_coef_cb store) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);

    if (mr_fvm_scalar_store_coef_cb_is_null(store)) {
        return MR_SUCCESS;
    }

    // TODO: Handle cross-derivative diffusion terms for non-orthogonal coordinates
    mr_float mul = 0.0f;
    mr_octree_cell_neighbor cell_neighbor = mr_octree_find_face_neighbor_cells(forest, cell_idx, dir);
    if (cell_neighbor.node_idx == MR_INVALID_INDEX || !mr_ocforest_is_node_active(forest, cell_neighbor.node_idx)) {
        mul = calc_transverse_der_mul(forest, cell_idx, dir);
        return interpolate_transverse_derivative(forest, cell_idx, dir, mul, store);
    }
    
    switch (cell_neighbor.type) {
        case MR_OCTREE_CELL_NEIGHBOR_NONE:
            return MR_FAILURE;

        case MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE:
            mul = calc_transverse_der_mul(forest, cell_idx, dir);

            return calc_equal_size_transverse_derivative(
                forest,
                cell_idx,
                cell_neighbor.neighbor_idx,
                dir,
                mul,
                store
            );

        case MR_OCTREE_CELL_NEIGHBOR_COARSER:
            mul = calc_transverse_der_mul(forest, cell_idx, dir);

            return calc_fine_coarse_transverse_derivative(
                forest,
                cell_idx,
                cell_neighbor.neighbor_idx,
                dir,
                mul,
                store
            );

        case MR_OCTREE_CELL_NEIGHBOR_FINER: ;
            mr_direction reflected_dir = mr_direction_reflect(dir);

            for (size_t i = 0; i < 4; ++i) {
                mul = -calc_transverse_der_mul(forest, cell_neighbor.neighbor_indices[i], reflected_dir);

                if (calc_fine_coarse_transverse_derivative(
                    forest,
                    cell_neighbor.neighbor_indices[i],
                    cell_idx,
                    reflected_dir,
                    mul,
                    store
                ) != MR_SUCCESS) {
                    return MR_FAILURE;
                }
            }
            return MR_SUCCESS;
    }

    return MR_FAILURE;
}

int mr_fvm_scalar_calc_transient_term(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_float64 step,
    mr_fvm_scalar_store_coef_cb store_implicit,
    mr_fvm_scalar_store_coef_cb store_rhs
) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);

    mr_float64 volume = mr_cell_volume(forest, cell_idx);
    mr_float64 coef = 1.0 / step * volume;

    if (!mr_fvm_scalar_store_coef_cb_is_null(store_implicit)) {
        store_implicit.fn(forest, cell_idx, coef, store_implicit.userdata);
    }

    if (!mr_fvm_scalar_store_coef_cb_is_null(store_rhs)) {
        store_rhs.fn(forest, cell_idx, coef, store_rhs.userdata);
    }

    return MR_SUCCESS;
}