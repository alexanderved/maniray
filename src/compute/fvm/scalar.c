#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "maniray/compute/fvm/scalar.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/fvm/cell.h"

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
    mr_fvm_scalar_store_coef_cb store_source
) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);
    assert(bc);

    switch (mr_boundary_condition_get_type(bc, cell_idx, dir)) {
        // Direchlet BC is implemented using second-order one-sided finite differences
        case MR_BC_DIRICHLET: ;
            mr_axis axis = mr_direction_get_axis(dir);
            mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

            // TODO: Move to helper function in cell.h
            mr_float face_center[] = { cell->x, cell->y, cell->z };
            face_center[axis] += mr_direction_get_sign_mul(dir) * cell->dim / 2.0f;

            // TODO: Create helper function to calculate derivatives
            mr_float sqrt_inv_coef = sqrtf(mr_manifold_inv_metric(forest->manifold, cell->chart_idx, face_center, axis, axis));
            mr_float area = mr_cell_face_area(forest, cell_idx, dir);

            if (!mr_fvm_scalar_store_coef_cb_is_null(store_implicit)) {
                mr_octree_cell_neighbor cell_neighbor = mr_octree_find_face_neighbor_cells(forest, cell_idx, mr_direction_reflect(dir));
                assert(cell_neighbor.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);

                // TODO: Handle non-orthogonal curvilinear coordinates
                mr_float cell_coef = -sqrt_inv_coef * 3.0f / cell->dim * area;
                mr_float neighbor_coef = sqrt_inv_coef * 1.0f / (3.0f * cell->dim) * area;

                store_implicit.fn(forest, cell_idx, cell_coef, store_implicit.userdata);
                store_implicit.fn(forest, cell_neighbor.neighbor_idx, neighbor_coef, store_implicit.userdata);
            }

            if (!mr_fvm_scalar_store_coef_cb_is_null(store_source)) {
                mr_float value = 0.0f;
                mr_boundary_condition_get_value(bc, cell_idx, dir, &value);

                mr_float coef = sqrt_inv_coef * 8.0f * value / (3.0f * cell->dim) * area;
                // Add minus sign because this term is moved over the equal sign to the rhs
                store_source.fn(forest, cell_idx, -coef, store_source.userdata);
            }

            break;

        case MR_BC_NEUMANN:
            if (!mr_fvm_scalar_store_coef_cb_is_null(store_source)) {
                mr_float area = mr_cell_face_area(forest, cell_idx, dir);
                mr_float value = 0.0f;
                mr_boundary_condition_get_value(bc, cell_idx, dir, &value);

                mr_float coef = value * area;
                // Add minus sign because this term is moved over the equal sign to the rhs
                store_source.fn(forest, cell_idx, -coef, store_source.userdata);
            }

            break;

        default:
            assert(false);
    }

    return MR_SUCCESS;
}

static int calc_equal_size_flux(
    mr_ocforest *forest,
    mr_int cell_idx,
    mr_int neighbor_cell_idx,
    mr_direction dir,
    mr_fvm_scalar_store_coef_cb store
) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_cell *neighbor_cell = mr_ocforest_get_cell(forest, neighbor_cell_idx);

    mr_axis axis = mr_direction_get_axis(dir);

    mr_float middle[] = {
        (cell->x + neighbor_cell->x) / 2.0f,
        (cell->y + neighbor_cell->y) / 2.0f,
        (cell->z + neighbor_cell->z) / 2.0f,
    };
    mr_float sqrt_inv_coef = sqrt(mr_manifold_inv_metric(forest->manifold, cell->chart_idx, middle, axis, axis));
    mr_float area = mr_cell_face_area(forest, cell_idx, dir);

    mr_float coef = sqrt_inv_coef / cell->dim * area;

    // TODO: Handle cross-derivative diffusion terms for non-orthogonal coordinates
    if (store.fn(forest, cell_idx, -coef, store.userdata) != MR_SUCCESS) {
        return MR_FAILURE;
    }
    return store.fn(forest, neighbor_cell_idx, coef, store.userdata);
}

typedef enum curr_cell_type {
    CURR_CELL_COARSE,
    CURR_CELL_FINE,
} curr_cell_type;

typedef struct ghost_cell_userdata {
    mr_fvm_scalar_store_coef_cb store;
    mr_float mul;
} ghost_cell_userdata;

static int store_ghost_cell_coefs(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    ghost_cell_userdata *gc_ud = userdata;
    return gc_ud->store.fn(forest, cell_idx, coef * gc_ud->mul, gc_ud->store.userdata);
}

static int calc_coarse_fine_flux(
    mr_ocforest *forest,
    mr_int coarse_cell_idx,
    mr_int fine_cell_idx,
    mr_direction dir,
    curr_cell_type type,
    mr_fvm_scalar_store_coef_cb store
) {
    mr_octree_cell *fine_cell = mr_ocforest_get_cell(forest, fine_cell_idx);

    mr_axis axis = mr_direction_get_axis(dir);
    mr_float middle[] = { fine_cell->x, fine_cell->y, fine_cell->z };
    middle[axis] += mr_direction_get_sign_mul(dir) * fine_cell->dim / 2.0f;

    mr_float sqrt_inv_coef = sqrt(mr_manifold_inv_metric(forest->manifold, fine_cell->chart_idx, middle, axis, axis));
    mr_float area = mr_cell_face_area(forest, fine_cell_idx, dir);

    mr_float coef = sqrt_inv_coef / fine_cell->dim * area * (type == CURR_CELL_FINE ? 1.0 : -1.0);
    if (store.fn(forest, fine_cell_idx, -coef, store.userdata) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    ghost_cell_userdata ud = { store, coef };
    return mr_fvm_calculate_ghost_cell(
        forest,
        coarse_cell_idx,
        fine_cell_idx,
        mr_fvm_interpolation_cb_create(store_ghost_cell_coefs, &ud)
    );
}

int mr_fvm_scalar_calc_internal_flux(mr_ocforest *forest, mr_int cell_idx, mr_direction dir, mr_fvm_scalar_store_coef_cb store) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);

    if (mr_fvm_scalar_store_coef_cb_is_null(store)) {
        return MR_SUCCESS;
    }

    mr_octree_cell_neighbor cell_neighbor = mr_octree_find_face_neighbor_cells(forest, cell_idx, dir);
    switch (cell_neighbor.type) {
        case MR_OCTREE_CELL_NEIGHBOR_NONE:
            return MR_FAILURE;
        
        case MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE:
            return calc_equal_size_flux(forest, cell_idx, cell_neighbor.neighbor_idx, dir, store);

        case MR_OCTREE_CELL_NEIGHBOR_COARSER:
            return calc_coarse_fine_flux(forest, cell_neighbor.neighbor_idx, cell_idx, dir, CURR_CELL_FINE, store);

        case MR_OCTREE_CELL_NEIGHBOR_FINER:
            for (size_t i = 0; i < 4; ++i) {
                if (calc_coarse_fine_flux(
                    forest,
                    cell_idx,
                    cell_neighbor.neighbor_indices[i],
                    mr_direction_reflect(dir),
                    CURR_CELL_COARSE,
                    store
                ) != MR_SUCCESS) {
                    return MR_FAILURE;
                }
            }
            return MR_SUCCESS;
    }

    return MR_FAILURE;
}