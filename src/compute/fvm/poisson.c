#include "maniray/utils/xmalloc.h"
#include "maniray/compute/fvm/grid.h"
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

typedef struct linear_system_data {
    mr_sparse_row *temp_row;
    mr_sparse_matrix_builder *matrix_builder;

    mr_float *source_terms;
} linear_system_data;

static void linear_system_data_create(mr_fvm_poisson *poisson, linear_system_data *data) {
    data->temp_row = mr_sparse_row_create();
    data->matrix_builder = mr_sparse_matrix_builder_create(poisson->code_map->len);

    data->source_terms = xcalloc(poisson->code_map->len, sizeof(mr_float));
}

static void linear_system_data_destroy(linear_system_data *data) {
    free(data->source_terms);
    mr_sparse_matrix_builder_destroy(data->matrix_builder);
    mr_sparse_row_destroy(data->temp_row);
}

static int fill_linear_system(mr_ocforest *forest, mr_int idx, void *userdata) {


    return MR_SUCCESS;
}

int mr_fvm_poisson_solve(mr_fvm_poisson *poisson) {
    if (!poisson) {
        return MR_FAILURE;
    }

    if (!poisson->code_map) {
        poisson->code_map = mr_code_map_create_from_ocforest(poisson->forest);
    }

    linear_system_data sys_data;
    linear_system_data_create(poisson, &sys_data);

    for (mr_index octree_idx = 0; octree_idx < poisson->forest->nb_roots; ++octree_idx) {
        int res = mr_octree_leaves_apply(
            poisson->forest,
            octree_idx,
            mr_octree_cond_cb_null(),
            mr_octree_apply_cb_create(fill_linear_system, &sys_data),
            false
        );

        if (res != MR_SUCCESS) {
            linear_system_data_destroy(&sys_data);

            return MR_FAILURE;
        }
    }

    linear_system_data_destroy(&sys_data);
}