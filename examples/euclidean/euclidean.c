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
#include "maniray/compute/fvm/general.h"
#include "maniray/compute/fvm/interpolation.h"

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

static int setup_boundary(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    for (mr_direction dir = MR_DIRECTION_MI_X; dir <= MR_DIRECTION_PL_Z; ++dir) {
        mr_int nidx = mr_octree_find_face_neighbor(forest, idx, dir);
        mr_octree_node *n = mr_ocforest_get_node(forest, nidx);
        if (!n) {
            discr_data->type = MR_CELL_TYPE_BOUNDARY;

            break;
        }
    }

    return MR_SUCCESS;
}

static bool point_refine(mr_ocforest *forest, mr_int idx, void *userdata) {
    mr_float *p = userdata;
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    return mr_norm_inf(p[0] - node->x, p[1] - node->y, p[2] - node->z) <= node->dim / 2.0f;
}

static bool area_refine(mr_ocforest *forest, mr_int idx, void *userdata) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    return node->y <= -1.0f; // && node->z >= 0.0;
}

static bool area_refine2(mr_ocforest *forest, mr_int idx, void *userdata) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    return node->y <= -1.0f; // && node->z >= 0.0;
}

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
    mr_ocforest *forest = mr_ocforest_create(manifold, descs, NB_ROOTS, (size_t[1]) { sizeof(mr_discretization_data) }, 1);

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);


    for (size_t i = 0; i < 1; ++i) {
        mr_octree_refine_all(forest, 0);
    }

    mr_float p[3] = { 0.5f, 0.5f, 0.5f };
    // mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_octree_cond_cb_create(point_refine, p), true);
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_octree_cond_cb_create(area_refine, NULL), false);
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_octree_cond_cb_create(area_refine2, NULL), false);
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_octree_cond_cb_create(area_refine2, NULL), false);
    mr_octree_balance(forest, 0);

    mr_int point_node_idx = mr_octree_locate_point(forest, 0, p);
    mr_octree_node *point_node = mr_ocforest_get_node(forest, point_node_idx);
    printf("Point Node %d: %f   (%f, %f, %f)\n",
        point_node_idx,
        point_node->dim,
        point_node->x,
        point_node->y,
        point_node->z
    );

    mr_direction vertex[] = { MR_DIRECTION_MI_X, MR_DIRECTION_MI_Y, MR_DIRECTION_PL_Z };
    mr_int neighbor_node_idx = mr_octree_find_vertex_neighbor(forest, point_node_idx, vertex);
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


    clock_gettime(CLOCK_MONOTONIC, &end);
    long long elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                           (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Refine + Balance: %.2f ms\n", (double)elapsed_us / 1000.0);


    clock_gettime(CLOCK_MONOTONIC, &start);

    mr_fvm_fit_grids_to_charts(forest);
    mr_octree_leaves_apply(forest, 0, mr_octree_cond_cb_null(), mr_octree_apply_cb_create(setup_boundary, NULL), false);
    mr_fvm_connect_overset_grids(forest);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                 (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Combine Grids: %.2f ms\n", (double)elapsed_us / 1000.0);


    printf("Ocforest size: %lu\n", mr_ocforest_size(forest));

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

    mr_isize octree_nodes_size = mr_ocforest_size(forest) * sizeof(mr_octree_node);
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