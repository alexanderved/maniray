#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "maniray/maniray.h"
#include "maniray/display/engine.h"
#include "maniray/display/camera.h"
#include "maniray/display/uniform_buffer.h"
#include "maniray/display/storage_buffer.h"
#include "maniray/utils/misc.h"
#include "maniray/compute/math.h"
#include "maniray/compute/manifold.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/fvm/grid.h"
#include "maniray/compute/fvm/cell.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/compute/fvm/poisson.h"

static bool chart_0_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);
    MR_UNUSED(p);

    return true;
}

static mr_manifold *setup_manifold() {
#define NB_CHARTS 1
    mr_chart_desc charts[NB_CHARTS] = {
        { .bounds = chart_0_bounds },
    };

    mr_transition_desc transitions[NB_CHARTS * NB_CHARTS] = {
        mr_transition_desc_create_self(),
    };

    return mr_manifold_create(3, NB_CHARTS, charts, transitions);
}

static int setup_boundary(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_discretization_data *discr_data = mr_ocforest_get_cell_extra(forest, cell_idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
        if (mr_is_boundary_cell(forest, cell_idx, dir)) {
            discr_data->type = MR_CELL_TYPE_BOUNDARY;

            break;
        }
    }

    return MR_SUCCESS;
}

static bool point_refine(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    mr_float *p = userdata;
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

    return mr_norm_inf(p[0] - cell->x, p[1] - cell->y, p[2] - cell->z) <= cell->dim / 2.0f;
}

static bool area_refine(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    return /* cell->x >= -1.0 && cell->y <= 1.0f && */ /* cell->z >= -1.0 && */ cell->z <= 0.0;

    // return cell->y >= -2.0f && cell->y <= 2.0f && mr_norm2_2d(cell->x, cell->z) <= 2.5f;
}

static int interpolation_test(mr_ocforest *forest, mr_int cell_idx, mr_float coef, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    printf("Interpolation Cell %d (Coef: %f): %f   (%f, %f, %f)\n",
        cell_idx,
        coef,
        cell->dim,
        cell->x,
        cell->y,
        cell->z
    );

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
    // return 0.0f;

    if (cell_idx == point_cell_idx) {
        return 1.0f / mr_cell_volume(poisson->forest, cell_idx) - 1.0f / 512.0f;
    }

    return -1.0f / 512.0f;
}

static mr_float value_min = 1e6f;
static mr_float value_max = -1e6f;

static int calc_bounds(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

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

// TODO: Investigate weird behaviour on partially refined grid
mr_boundary_condition *setup_bc(mr_ocforest *forest) {
    mr_boundary_condition_type bc_type = MR_BC_DIRICHLET;
    mr_boundary_condition_fn bc_fn = bc_zero;

    mr_boundary_condition_type bc_types[][MR_NB_DIRECTIONS] = { { bc_type, bc_type, bc_type, bc_type, bc_type, bc_type } };
    mr_boundary_condition_fn bc_fns[] = { bc_fn, bc_fn, bc_fn, bc_fn, bc_fn, bc_fn };
    mr_boundary_condition *bc = mr_boundary_condition_create(forest, bc_types, bc_fns);

    return bc;
}

#define INIT_TIMER() struct timespec __start, __end;
#define START_TIMER() clock_gettime(CLOCK_MONOTONIC, &__start);
#define STOP_TIMER(str) \
    do { \
        clock_gettime(CLOCK_MONOTONIC, &__end); \
        long long __elapsed_us = (__end.tv_sec - __start.tv_sec) * 1000000LL + \
                            (__end.tv_nsec - __start.tv_nsec) / 1000; \
        printf("%s: %.2f ms\n", (str), (double)__elapsed_us / 1000.0); \
    } while (0)

mr_ocforest *setup_ocforest(mr_manifold *manifold) {
#define NB_ROOTS 1
    mr_octree_root_desc descs[NB_ROOTS] = {
        {
            .flags = 0,
            .chart_idx = 0,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 8.0f,
        },
    };

    mr_fvm_poisson *poisson = mr_fvm_poisson_create();
    mr_ocforest *forest = mr_fvm_poisson_ocforest_initialize(poisson, manifold, descs, NB_ROOTS);


    INIT_TIMER();
    START_TIMER();

    mr_octree_refine_all(forest, 0, 3);

    mr_float p[3] = { 0.5f, 0.5f, -0.5f };
    mr_octree_refine(forest, 0, mr_octree_cond_cb_create(point_refine, p), false);
    // mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_octree_cond_cb_create(point_refine, (mr_float[]) { -0.5f, -0.5f, -0.5f }), false);
    mr_octree_refine(forest, 0, mr_octree_cond_cb_create(area_refine, NULL), false);
    mr_octree_balance(forest, 0);

    mr_octree_cells_apply(forest, 0, mr_octree_apply_cb_create(setup_boundary, NULL));

    STOP_TIMER("Refine + Balance");


    mr_fvm_poisson_ocforest_finalize(poisson);
    mr_fvm_poisson_set_boundary_condition(poisson, setup_bc(forest));
    mr_fvm_poisson_set_source_term_fn(poisson, source_test);

    point_cell_idx = mr_octree_locate_point_in_cell(forest, 0, (mr_float[]) { 0.1f, -0.1f, -0.5f });


    START_TIMER();

    mr_fvm_poisson_build_discretization_matrix(poisson);
    mr_fvm_poisson_build_source_terms(poisson);

    STOP_TIMER("Build Equation");


    START_TIMER();

    mr_fvm_poisson_solve(poisson);

    STOP_TIMER("Solve");


    mr_octree_cells_apply(forest, 0, mr_octree_apply_cb_create(calc_bounds, NULL));
    mr_octree_cells_apply(forest, 0, mr_octree_apply_cb_create(normalize_eq_res, NULL));

    printf("MIN / MAX: %f / %f\n", value_min, value_max);


    START_TIMER();

    // mr_fvm_fit_grids_to_charts(forest);
    // mr_fvm_connect_overset_grids(forest);

    STOP_TIMER("Combine Grids");


    printf("Number of nodes: %lu\n", mr_ocforest_nb_nodes_upper_bound(forest));
    printf("Number of cells: %lu\n", mr_ocforest_nb_cells_upper_bound(forest));

    return forest;
}

typedef struct update_userdata {
    mr_window *window;
    mr_camera *camera;
    mr_uniform_buffer *camera_buffer;
} update_userdata;

void update_func(mr_engine *engine, void *userdata) {
    MR_UNUSED(engine);

    update_userdata *data = userdata;

    mr_window *window = data->window;
    mr_camera *camera = data->camera;
    mr_uniform_buffer *camera_buffer = data->camera_buffer;

    mr_camera_info camera_info = mr_camera_get_info(camera);
    mr_uniform_buffer_fill(camera_buffer, 0, sizeof(camera_info), &camera_info);

    mr_camera_proccess_movement(camera, window);

    // printf("FPS: %f\n", 1.0 / mr_window_get_delta_time(window));
}

int run_display() {
    mr_window_init();

    mr_window *window = mr_window_create(MR_WINDOW_DEFAULT_WIDTH, MR_WINDOW_DEFAULT_HEIGHT, "ManiRay");
    if (!window) {
        printf("Failed to create a window\n");
        mr_window_terminate();

        return MR_FAILURE;
    }

    if (mr_engine_init(window) != MR_SUCCESS) {
        printf("Failed to initialize GLAD\n");
        mr_window_terminate();

        return MR_FAILURE;
    }
    
    mr_engine *engine = mr_engine_create(window);

    mr_camera *camera = mr_camera_create(window, (vec3){ 0.0, 0.0, -3.0 }, 90.0f, 1.0f);
    if (!camera) {
        mr_engine_destroy(engine);
        mr_window_destroy(window);
        mr_window_terminate();

        return MR_FAILURE;
    }

    mr_uniform_buffer *camera_buffer = mr_uniform_buffer_create(2, sizeof(mr_camera_info), MR_DYNAMIC_DRAW);
    if (!camera_buffer) {
        mr_camera_destroy(camera);
        mr_engine_destroy(engine);
        mr_window_destroy(window);
        mr_window_terminate();

        return MR_FAILURE;
    }

    mr_camera_info camera_info = mr_camera_get_info(camera);
    mr_uniform_buffer_fill(camera_buffer, 0, sizeof(camera_info), &camera_info);

    mr_manifold *manifold = setup_manifold();
    mr_ocforest *forest = setup_ocforest(manifold);

    // Octree nodes buffer
    mr_isize octree_nodes_size = mr_ocforest_nb_nodes_upper_bound(forest) * sizeof(mr_octree_node);
    mr_isize octree_buffer_size = sizeof(mr_uint) + octree_nodes_size;
    mr_storage_buffer *octree_buffer = mr_storage_buffer_create(0);
    mr_storage_buffer_alloc(octree_buffer, octree_buffer_size, MR_STATIC_DRAW, NULL);

    char *octree_buffer_ptr = mr_storage_buffer_map(octree_buffer, 0, octree_buffer_size, MR_STORAGE_BUFFER_FLAG_WRITE);
    mr_uint nb_roots = (mr_uint)forest->nb_roots;
    memcpy(octree_buffer_ptr, &nb_roots, sizeof(mr_uint));
    memcpy(octree_buffer_ptr + sizeof(mr_uint), mr_ocforest_get_node_array(forest), octree_nodes_size);
    mr_storage_buffer_unmap(octree_buffer);

    // PDE solution buffer
    mr_isize solution_size = mr_ocforest_nb_cells_upper_bound(forest) * sizeof(mr_fvm_poisson_solution);
    mr_storage_buffer *solution_buffer = mr_storage_buffer_create(1);
    mr_storage_buffer_alloc(solution_buffer, solution_size, MR_STATIC_DRAW, mr_ocforest_get_cell_extra_array(forest, MR_POISSON_SOLUTION_EXTRA_FIELD));

    update_userdata update_data = { window, camera, camera_buffer };
    
    mr_shader_source compute_shader_source = mr_shader_source_read("../examples/euclidean/euclidean_slice.comp", MR_COMPUTE);
    mr_program *compute_program = mr_program_create((mr_shader_source[]){ compute_shader_source }, 1);
    mr_shader_source_destroy(compute_shader_source);

    if (!compute_program) {
        mr_storage_buffer_destroy(octree_buffer);
        mr_ocforest_destroy(forest);

        mr_uniform_buffer_destroy(camera_buffer);
        mr_camera_destroy(camera);

        mr_engine_destroy(engine);
        mr_window_destroy(window);
        mr_window_terminate();

        return MR_FAILURE;
    }
    
    // mr_window_set_cursor_hidden(window, true);
    mr_engine_run(engine, compute_program, update_func, &update_data);

    mr_program_destroy(compute_program);

    mr_storage_buffer_destroy(solution_buffer);
    mr_storage_buffer_destroy(octree_buffer);
    mr_ocforest_destroy(forest);

    mr_uniform_buffer_destroy(camera_buffer);
    mr_camera_destroy(camera);

    mr_engine_destroy(engine);
    mr_window_destroy(window);
    mr_window_terminate();

    return MR_SUCCESS;
}

int main(int argc, char *argv[]) {
    mr_initialize(&argc, &argv);

#define DISPLAY
#ifdef DISPLAY
    run_display();
#else
    mr_manifold *manifold = setup_manifold();
    mr_ocforest *forest = setup_ocforest(manifold);
#endif

    mr_finalize();

    return 0;
}