#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "maniray/display/engine.h"
#include "maniray/display/camera.h"
#include "maniray/display/uniform_buffer.h"
#include "maniray/display/storage_buffer.h"
#include "maniray/utils/misc.h"
#include "maniray/compute/math.h"
#include "maniray/compute/manifold.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/fvm/grid.h"
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

static bool point_refine(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    mr_float *p = userdata;
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);

    printf("%f %f\n", mr_norm_inf(p[0] - cell->x, p[1] - cell->y, p[2] - cell->z), cell->dim / 2.0f);

    return mr_norm_inf(p[0] - cell->x, p[1] - cell->y, p[2] - cell->z) <= cell->dim / 2.0f;
}

static bool area_refine(mr_ocforest *forest, mr_int cell_idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    return cell->y <= -1.0f; // && cell->z >= 0.0;
}

#if 0
static int interpolation_test(mr_ocforest *forest, mr_int node_idx, mr_float coef, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);
    printf("Interpolation Node %d (Coef: %f): %f   (%f, %f, %f)\n",
        node_idx,
        coef,
        node->dim,
        node->x,
        node->y,
        node->z
    );

    return MR_SUCCESS;
}
#endif

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

    mr_fvm_poisson *poisson = mr_fvm_poisson_create(manifold, descs, NB_ROOTS, NULL);
    mr_ocforest *forest = poisson->forest;

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);


    mr_octree_refine_all(forest, 0, 2);

    mr_float p[3] = { 0.5f, 0.5f, -0.5f };
    mr_octree_refine(forest, 0, mr_octree_cond_cb_create(point_refine, p), false);
    /* mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_octree_cond_cb_create(point_refine, (mr_float[]) { -0.5f, -0.5f, -0.5f }), false); */
    /* mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_octree_cond_cb_create(area_refine, NULL), false); */
    mr_octree_balance(forest, 0);

    mr_octree_cells_apply(forest, 0, mr_octree_apply_cb_create(setup_boundary, NULL));
    // mr_fvm_poisson_build_discretization_matrix(poisson);

    clock_gettime(CLOCK_MONOTONIC, &end);
    long long elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                           (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Refine + Balance: %.2f ms\n", (double)elapsed_us / 1000.0);


    mr_int point_node_idx = mr_octree_locate_point_in_leaf(forest, 0, (mr_float[]) { 0.5f, 0.5f, -0.5f });
    mr_octree_node *point_node = mr_ocforest_get_node(forest, point_node_idx);
    printf("Point Node %d: %f   (%f, %f, %f)\n",
        point_node_idx,
        point_node->dim,
        point_node->x,
        point_node->y,
        point_node->z
    );

#if 0
    mr_direction edge[] = { MR_DIRECTION_MI_X, MR_DIRECTION_MI_Y };
    mr_int neighbor_node_idx = mr_octree_find_edge_neighbor(forest, point_node_idx, edge);
    mr_octree_node *neighbor_node = mr_ocforest_get_node(forest, neighbor_node_idx);
    if (neighbor_node) {
        printf("Neighbor Node %d: %f   (%f, %f, %f)\n",
            neighbor_node_idx,
            neighbor_node->dim,
            neighbor_node->x,
            neighbor_node->y,
            neighbor_node->z
        );
    }

    mr_fvm_calculate_ghost_cell(
        forest,
        point_node_idx,
        (mr_direction[]) { MR_DIRECTION_MI_X, MR_DIRECTION_MI_Y, MR_DIRECTION_MI_Z },
        mr_fvm_interpolation_cb_create(interpolation_test, NULL)
    );
#endif


    clock_gettime(CLOCK_MONOTONIC, &start);

    // mr_fvm_fit_grids_to_charts(forest);
    // mr_fvm_connect_overset_grids(forest);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                 (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Combine Grids: %.2f ms\n", (double)elapsed_us / 1000.0);


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

        return -1;
    }

    if (!mr_engine_init(window)) {
        printf("Failed to initialize GLAD\n");
        mr_window_terminate();

        return -1;
    }
    
    mr_engine *engine = mr_engine_create(window);

    mr_camera *camera = mr_camera_create(window, (vec3){ 0.0, 0.0, -3.0 }, 90.0f, 1.0f);
    if (!camera) {
        mr_engine_destroy(engine);
        mr_window_destroy(window);
        mr_window_terminate();

        return -1;
    }

    mr_uniform_buffer *camera_buffer = mr_uniform_buffer_create(2, sizeof(mr_camera_info), MR_DYNAMIC_DRAW);
    if (!camera_buffer) {
        mr_camera_destroy(camera);
        mr_engine_destroy(engine);
        mr_window_destroy(window);
        mr_window_terminate();

        return -1;
    }

    mr_camera_info camera_info = mr_camera_get_info(camera);
    mr_uniform_buffer_fill(camera_buffer, 0, sizeof(camera_info), &camera_info);

    mr_manifold *manifold = setup_manifold();
    mr_ocforest *forest = setup_ocforest(manifold);

    mr_isize octree_nodes_size = mr_ocforest_nb_nodes_upper_bound(forest) * sizeof(mr_octree_node);
    mr_isize octree_buffer_size = sizeof(mr_uint) + octree_nodes_size;
    mr_storage_buffer *octree_buffer = mr_storage_buffer_create(0);
    mr_storage_buffer_alloc(octree_buffer, octree_buffer_size, MR_STATIC_DRAW, NULL);

    char *octree_buffer_ptr = mr_storage_buffer_map(octree_buffer, 0, octree_buffer_size, MR_STORAGE_BUFFER_FLAG_WRITE);
    mr_uint nb_roots = (mr_uint)forest->nb_roots;
    memcpy(octree_buffer_ptr, &nb_roots, sizeof(mr_uint));
    memcpy(octree_buffer_ptr + sizeof(mr_uint), mr_ocforest_get_node_array(forest), octree_nodes_size);
    mr_storage_buffer_unmap(octree_buffer);

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

        return -1;
    }
    
    // mr_window_set_cursor_hidden(window, true);
    mr_engine_run(engine, compute_program, update_func, &update_data);

    mr_program_destroy(compute_program);

    mr_storage_buffer_destroy(octree_buffer);
    mr_ocforest_destroy(forest);

    mr_uniform_buffer_destroy(camera_buffer);
    mr_camera_destroy(camera);

    mr_engine_destroy(engine);
    mr_window_destroy(window);
    mr_window_terminate();

    return 0;
}

int main() {
    run_display();

    return 0;
}