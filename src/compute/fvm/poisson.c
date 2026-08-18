#include <stdio.h>
#include <tgmath.h>

#include "maniray/utils/xmalloc.h"
#include "maniray/compute/fvm/cell.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/compute/fvm/poisson.h"
#include "maniray/compute/sparse_matrix.h"

mr_fvm_poisson *mr_fvm_poisson_create(
    mr_manifold *manifold,
    mr_octree_root_desc roots[],
    size_t nb_roots,
    mr_fvm_poisson_source_fn source_fn
) {
    mr_fvm_poisson *poisson = xmalloc(sizeof(mr_fvm_poisson));

    poisson->forest = mr_ocforest_create(
        manifold,
        roots,
        nb_roots,
        (size_t[]) { sizeof(mr_discretization_data), sizeof(mr_fvm_poisson_solution) },
        2
    );
    poisson->code_map = NULL;

    poisson->source_fn = source_fn;

    return poisson;
}

void mr_fvm_poisson_destroy(mr_fvm_poisson *poisson) {
    if (!poisson) {
        return;
    }

    mr_code_map_destroy(poisson->code_map);
    mr_ocforest_destroy(poisson->forest);

    free(poisson);
}

mr_ocforest *mr_fvm_poisson_get_ocforest(mr_fvm_poisson *poisson) {
    return poisson ? poisson->forest : NULL;
}

// TODO: Make only for matrix (exclude source terms). Rename to fit new purpose
typedef struct linear_system_data {
    mr_fvm_poisson *poisson;

    mr_sparse_row *temp_row;
    mr_sparse_matrix_builder *matrix_builder;

    mr_float *source_terms;
} linear_system_data;

static void linear_system_data_create(mr_fvm_poisson *poisson, linear_system_data *data) {
    data->poisson = poisson;

    data->temp_row = mr_sparse_row_create();
    data->matrix_builder = mr_sparse_matrix_builder_create(poisson->code_map->len);

    data->source_terms = xcalloc(poisson->code_map->len, sizeof(mr_float));
}

static void linear_system_data_destroy(linear_system_data *data) {
    free(data->source_terms);
    mr_sparse_matrix_builder_destroy(data->matrix_builder);
    mr_sparse_row_destroy(data->temp_row);
}

static void add_to_cell_matrix_coef(mr_ocforest *forest, mr_int idx, mr_float value, linear_system_data *sys_data) {
    mr_int code = mr_ocforest_get_code(forest, idx);
    size_t col = mr_code_map_get_index(sys_data->poisson->code_map, code);

    mr_float prev = mr_sparse_row_get(sys_data->temp_row, col);
    mr_sparse_row_set(sys_data->temp_row, col, prev + value);
}

static int fill_inactive_cell_row(mr_ocforest *forest, mr_int idx, linear_system_data *sys_data) {
    add_to_cell_matrix_coef(forest, idx, 1.0f, sys_data);

    return MR_SUCCESS;
}

static int get_interp_point(mr_ocforest *forest, mr_int idx, mr_float p[MR_NB_AXES]) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);

    mr_int donor_root_node_idx = forest->roots[discr_data->donor_root_idx].node_idx;
    mr_uint donor_chart = mr_ocforest_get_node(forest, donor_root_node_idx)->chart_idx;

    mr_float center[MR_NB_AXES] = { node->x, node->y, node->z };
    return mr_manifold_transition(forest->manifold, node->chart_idx, donor_chart, p, center);
}

static int write_interp_coef(mr_ocforest *forest, mr_int idx, mr_float coef, void *userdata) {
    linear_system_data *sys_data = userdata;
    add_to_cell_matrix_coef(forest, idx, -coef, sys_data);

    return MR_SUCCESS;
}

static int fill_interp_cell_row(mr_ocforest *forest, mr_int idx, linear_system_data *sys_data) {
    add_to_cell_matrix_coef(forest, idx, 1.0f, sys_data);

    mr_float p[MR_NB_AXES] = { 0.0f };
    if (get_interp_point(forest, idx, p) != MR_SUCCESS) {
        return MR_FAILURE;
    }

    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);
    return mr_fvm_perform_interpolation(
        forest,
        discr_data->donor_root_idx,
        p,
        mr_fvm_interpolation_cb_create(write_interp_coef, sys_data)
    );
}

static int calc_boundary_flux_implicit_terms(mr_ocforest *forest, mr_int idx, mr_direction dir, linear_system_data *sys_data) {
    MR_UNUSED(forest);
    MR_UNUSED(idx);
    MR_UNUSED(dir);
    MR_UNUSED(sys_data);

    // TODO: Implement different boundary conditions
    // Currently using Neumann boundary conditions with flux = 0

    return MR_SUCCESS;
}

static int calc_same_size_flux(mr_ocforest *forest, mr_int idx, mr_int nidx, mr_direction dir, linear_system_data *sys_data) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_octree_node *neighbor_node = mr_ocforest_get_node(forest, nidx);

    mr_axis axis = mr_direction_get_axis(dir);

    mr_float middle[] = {
        (node->x + neighbor_node->x) / 2.0f,
        (node->y + neighbor_node->y) / 2.0f,
        (node->z + neighbor_node->z) / 2.0f,
    };
    mr_float sqrt_inv_coef = sqrt(mr_manifold_inv_metric(forest->manifold, node->chart_idx, middle, axis, axis));
    mr_float area = mr_cell_face_area(forest, idx, dir);

    mr_float coef = sqrt_inv_coef / node->dim * area;

    // TODO: Handle cross-derivative diffusion terms for non-orthogonal coordinates
    add_to_cell_matrix_coef(forest, idx, -coef, sys_data);
    add_to_cell_matrix_coef(forest, nidx, coef, sys_data);

    return MR_SUCCESS;
}

typedef struct ghost_cell_userdata {
    linear_system_data *sys_data;
    mr_float mul;
} ghost_cell_userdata;

static int write_ghost_cell_coef(mr_ocforest *forest, mr_int idx, mr_float coef, void *userdata) {
    ghost_cell_userdata *gc_ud = userdata;
    add_to_cell_matrix_coef(forest, idx, coef * gc_ud->mul, gc_ud->sys_data);

    return MR_SUCCESS;
}

// TODO(?): Combine with calc_same_size_flux
static int calc_coarse_fine_flux(mr_ocforest *forest, mr_int idx, mr_int nidx, mr_direction dir, mr_sign sign, linear_system_data *sys_data) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_octree_node *parent = mr_ocforest_get_node(forest, node->parent);

    mr_axis axis = mr_direction_get_axis(dir);
    mr_float middle[] = { node->x, node->y, node->z };
    middle[axis] += mr_direction_get_sign_mul(dir) * node->dim / 2.0f;

    mr_float sqrt_inv_coef = sqrt(mr_manifold_inv_metric(forest->manifold, node->chart_idx, middle, axis, axis));
    mr_float area = mr_cell_face_area(forest, idx, dir);
    mr_float coef = sqrt_inv_coef / node->dim * area * mr_sign_to_mul(sign);

    add_to_cell_matrix_coef(forest, idx, -coef, sys_data);

    mr_int local_idx = idx - parent->first_child;
    mr_direction ghost_cell_dir[MR_ADJACENCY_VERTEX] = { 0 };
    mr_local_idx_to_direction(local_idx, ghost_cell_dir);

    ghost_cell_userdata ud = { sys_data, coef };
    return mr_fvm_calculate_ghost_cell(
        forest,
        nidx,
        ghost_cell_dir,
        mr_fvm_interpolation_cb_create(write_ghost_cell_coef, &ud)
    );
}

static int fill_discr_cell_row(mr_ocforest *forest, mr_int idx, linear_system_data *sys_data) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);

    for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
        mr_int nidx = mr_octree_find_face_neighbor(forest, idx, dir);
        if (nidx == MR_INVALID_INDEX) {
            if (discr_data->type != MR_CELL_TYPE_BOUNDARY) {
                return MR_FAILURE;
            }

            if (calc_boundary_flux_implicit_terms(forest, idx, dir, sys_data) != MR_SUCCESS) {
                return MR_FAILURE;
            }

            continue;
        }

        mr_octree_node *neighbor_node = mr_ocforest_get_node(forest, nidx);
        if (neighbor_node->flags & MR_OCTREE_NODE_FLAG_LEAF && node->level == neighbor_node->level) {
            if (calc_same_size_flux(forest, idx, nidx, dir, sys_data) != MR_SUCCESS) {
                return MR_FAILURE;
            }
        } else if (neighbor_node->flags & MR_OCTREE_NODE_FLAG_LEAF) {
            if (calc_coarse_fine_flux(forest, idx, nidx, dir, MR_SIGN_PLUS, sys_data) != MR_SUCCESS) {
                return MR_FAILURE;
            }
        } else {
            mr_direction reflected_dir = mr_direction_reflect(dir);
            for (mr_int local_idx = 0; local_idx < MR_OCTREE_NB_CHILDREN; ++local_idx) {
                if (!mr_is_local_idx_face_adjacent(local_idx, reflected_dir)) {
                    continue;
                }

                mr_int child_idx = neighbor_node->first_child + local_idx;
                if (calc_coarse_fine_flux(forest, child_idx, idx, reflected_dir, MR_SIGN_MINUS, sys_data) != MR_SUCCESS) {
                    return MR_FAILURE;
                }
            }
        }
    }

    return MR_SUCCESS;
}

static int fill_discr_matrix(mr_ocforest *forest, mr_int idx, void *userdata) {
    linear_system_data *sys_data = userdata;

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);

    int res = MR_SUCCESS;
    if (!(node->flags & MR_OCTREE_NODE_FLAG_ACTIVE)) {
        res = fill_inactive_cell_row(forest, idx, sys_data);
    } else if (discr_data->type == MR_CELL_TYPE_INTERPOLATION) {
        res = fill_interp_cell_row(forest, idx, sys_data);
    } else {
        res = fill_discr_cell_row(forest, idx, sys_data);
    }

    if (res != MR_SUCCESS) {
        return res;
    }

    mr_sparse_matrix_builder_add_row(sys_data->matrix_builder, sys_data->temp_row);
    mr_sparse_row_clear(sys_data->temp_row);

    return MR_SUCCESS;
}

int mr_fvm_poisson_build_discretization_matrix(mr_fvm_poisson *poisson) {
    if (!poisson) {
        return MR_FAILURE;
    }

    if (!poisson->code_map) {
        poisson->code_map = mr_code_map_create_from_ocforest(poisson->forest);
    }

    linear_system_data sys_data;
    linear_system_data_create(poisson, &sys_data);

    for (mr_index octree_idx = 0; (size_t)octree_idx < poisson->forest->nb_roots; ++octree_idx) {
        int res = mr_octree_leaves_apply(
            poisson->forest,
            octree_idx,
            mr_octree_cond_cb_null(),
            mr_octree_apply_cb_create(fill_discr_matrix, &sys_data),
            false
        );

        if (res != MR_SUCCESS) {
            linear_system_data_destroy(&sys_data);

            return MR_FAILURE;
        }
    }

#if 0
    for (size_t i = 0; i < sys_data.matrix_builder->dim; ++i) {
        printf("Row %lu: ", i);
        for (size_t j = sys_data.matrix_builder->rows[i]; j < sys_data.matrix_builder->rows[i + 1]; ++j) {
            printf("%lu: %f\t", sys_data.matrix_builder->cols[j], sys_data.matrix_builder->values[j]);
        }
        printf("\n");
    }
#endif

    linear_system_data_destroy(&sys_data);

    return MR_SUCCESS;
}