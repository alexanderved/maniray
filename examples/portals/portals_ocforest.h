#include <time.h>

#include "maniray/utils/misc.h"
#include "maniray/utils/types.h"
#include "maniray/compute/math.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/fvm/poisson.h"
#include "maniray/compute/fvm/cell.h"

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

static mr_int point_cell_idx = 0;

static void bc_zero(mr_boundary_condition *bc, mr_int cell_idx, mr_direction dir, void *out) {
    MR_UNUSED(bc);
    MR_UNUSED(cell_idx);
    MR_UNUSED(dir);

    *(mr_float *)out = 0.0;
}

static mr_float source_test(mr_fvm_poisson *poisson, mr_int cell_idx) {
    mr_float density = 8.0f;

    if (cell_idx == point_cell_idx) {
        return density / mr_cell_volume(poisson->forest, cell_idx);
    }

    return 0.0;
}

static mr_float value_min = 1e6f;
static mr_float value_max = -1e6f;

static int calc_bounds(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type == MR_CELL_TYPE_EXTERIOR) {
        return MR_SUCCESS;
    }

    mr_fvm_poisson_solution *sol = mr_ocforest_get_cell_extra(forest, cell_idx, MR_POISSON_SOLUTION_EXTRA_FIELD);
    value_min = MR_MIN(value_min, sol->value);
    value_max = MR_MAX(value_max, sol->value);

    return MR_SUCCESS;
}

static int normalize_eq_res(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_fvm_poisson_solution *sol = mr_ocforest_get_cell_extra(forest, cell_idx, MR_POISSON_SOLUTION_EXTRA_FIELD);

    if (MR_ABS(value_max - value_min) > 1.0e-3f) {
        sol->value = (sol->value - value_min) / (value_max - value_min);
    } else {
        sol->value = 1.0f;
    }

    return MR_SUCCESS;
}

mr_boundary_condition *setup_bc(mr_ocforest *forest) {
    mr_boundary_condition_type bc_types[][MR_NB_DIRECTIONS] = {
        { MR_BC_DIRICHLET, MR_BC_DIRICHLET, MR_BC_DIRICHLET, MR_BC_DIRICHLET, MR_BC_DIRICHLET, MR_BC_DIRICHLET },
        { MR_BC_NONE, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE },
        { MR_BC_NONE, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE },
        { MR_BC_DIRICHLET, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE, MR_BC_NONE },
    };
    mr_boundary_condition_fn bc_fns[] = { bc_zero, bc_zero, bc_zero, bc_zero };
    mr_boundary_condition *bc = mr_boundary_condition_create(forest, bc_types, bc_fns);

    return bc;
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
    
    mr_fvm_poisson *poisson = mr_fvm_poisson_create();
    mr_ocforest *forest = mr_fvm_poisson_ocforest_initialize(poisson, manifold, descs, NB_ROOTS);


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





    mr_fvm_poisson_ocforest_finalize(poisson);
    mr_fvm_poisson_set_boundary_condition(poisson, setup_bc(forest));
    mr_fvm_poisson_set_source_term_fn(poisson, source_test);

    point_cell_idx = mr_octree_locate_point_in_cell(forest, 0, (mr_float[]) { 0.0f, 0.0f, 0.0f });



    mr_fvm_poisson_build_discretization_matrix(poisson);
    mr_fvm_poisson_build_source_terms(poisson);

    mr_linear_system_solver_set_options(poisson->solver, MR_SOLVER_BICGSTAB, MR_PRECON_SSOR, 1.0e-8);
    mr_linear_system_solver_print_debug_info(poisson->solver);
    mr_fvm_poisson_solve(poisson);


    for (size_t i = 0; i < forest->nb_roots; ++i) {
        mr_octree_cells_apply(forest, i, mr_octree_apply_cb_create(calc_bounds, NULL));
    }

    for (size_t i = 0; i < forest->nb_roots; ++i) {
        mr_octree_cells_apply(forest, i, mr_octree_apply_cb_create(normalize_eq_res, NULL));
    }


    printf("MIN / MAX: %f / %f\n", value_min, value_max);


    printf("Ocforest number of nodes: %lu\n", mr_ocforest_nb_nodes_upper_bound(forest));
    printf("Ocforest number of cells: %lu\n", mr_ocforest_nb_cells_upper_bound(forest));

    return forest;
}