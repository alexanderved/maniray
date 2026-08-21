#include <time.h>

#include "maniray/utils/misc.h"
#include "maniray/utils/types.h"
#include "maniray/compute/math.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/fvm/grid.h"

static bool adaptive_refine(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_float p[3] = { cell->x, cell->y, cell->z };

    return p[1] < 1.25f && p[1] > -1.25f && (mr_norm2_2d(p[0] + 2.0f, p[2]) < 1.8f || mr_norm2_2d(p[0] - 2.0f, p[2]) < 1.8f);
}

static int setup_boundary(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    mr_int local_idx = cell_idx - node->first_child;
    for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
        if (!mr_is_cell_local_idx_face_adjacent(local_idx, dir)) {
            continue;
        }

        mr_int nidx = mr_octree_find_face_neighbor_node(forest, cell->parent, dir);
        if (nidx == MR_INVALID_INDEX) {
            discr_data->type = MR_CELL_TYPE_BOUNDARY;

            break;
        }
    }

    return MR_SUCCESS;
}

static int setup_boundary_chart3(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    mr_direction dir = MR_DIRECTION_MI_X;
    mr_int local_idx = cell_idx - node->first_child;

    if (!mr_is_cell_local_idx_face_adjacent(local_idx, dir)) {
        return MR_SUCCESS;
    }

    mr_int nidx = mr_octree_find_face_neighbor_node(forest, cell->parent, dir);
    if (nidx == MR_INVALID_INDEX) {
        discr_data->type = MR_CELL_TYPE_BOUNDARY;
    }

    return MR_SUCCESS;
}

mr_ocforest *setup_ocforest(mr_manifold *manifold) {
#define NB_ROOTS 4
    mr_octree_root_desc descs[NB_ROOTS] = {
        {
            .flags = 0,
            .chart_idx = 0,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 8.0f,
        },
        {
            .flags = 0,
            .chart_idx = 1,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 2.0f,
        },
        {
            .flags = 0,
            .chart_idx = 2,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 2.0f,
        },
        {
            .flags = MR_OCTREE_FLAG_PERIODIC_Y | MR_OCTREE_FLAG_PERIODIC_Z,
            .chart_idx = 3,
            .x = 0.375f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 0.5f,
        },
    };
    mr_ocforest *forest = mr_ocforest_create(manifold, descs, NB_ROOTS, (size_t[1]) { sizeof(mr_discretization_data) }, 1);


    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    size_t min_refinement_level = MR_OCTREE_MAX_LEVEL - 2;
    mr_octree_refine_all(forest, 0, min_refinement_level);
    mr_octree_refine_all(forest, 1, min_refinement_level);
    mr_octree_refine_all(forest, 2, min_refinement_level);
    mr_octree_refine_all(forest, 3, min_refinement_level);

    mr_octree_refine(forest, 0, mr_octree_cond_cb_create(adaptive_refine, NULL), true);
    mr_octree_balance(forest, 0);

    clock_gettime(CLOCK_MONOTONIC, &end);
    long long elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                           (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Refine + Balance: %.2f ms\n", (double)elapsed_us / 1000.0);

    


    clock_gettime(CLOCK_MONOTONIC, &start);

    mr_fvm_fit_grids_to_charts(forest);
    mr_octree_cells_apply(forest, 0, mr_octree_apply_cb_create(setup_boundary, NULL));
    mr_octree_cells_apply(forest, 3, mr_octree_apply_cb_create(setup_boundary_chart3, NULL));
    mr_fvm_connect_overset_grids(forest);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                 (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Combine Grids: %.2f ms\n", (double)elapsed_us / 1000.0);


    printf("Ocforest number of nodes: %lu\n", mr_ocforest_nb_nodes_upper_bound(forest));
    printf("Ocforest number of cells: %lu\n", mr_ocforest_nb_cells_upper_bound(forest));

    return forest;
}