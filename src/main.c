#include <stdio.h>
#include <time.h>

#include "maniray/display/engine.h"
#include "maniray/display/camera.h"
#include "maniray/display/uniform_buffer.h"
#include "maniray/display/storage_buffer.h"
#include "maniray/utils/misc.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/sparse_matrix.h"

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

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

static float l2_2d(float x, float y) {
    return sqrtf(x * x + y * y);
}

static float l2(float x, float y, float z) {
    return sqrtf(x * x + y * y + z * z);
}

static float linf(float x, float y, float z) {
    return MAX(MAX(fabsf(x), fabsf(y)), fabsf(z));
}

static bool chart0_bounds(const mr_manifold *manifold, const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(manifold);
    MR_UNUSED(chart);

    return p[1] > 0.25f || p[1] < -0.25f || (l2_2d(p[0] + 2.0f, p[2]) > 1.0f && l2_2d(p[0] - 2.0f, p[2]) > 1.0f);
}

mr_manifold *setup_manifold() {
#define NB_CHARTS 1
    mr_chart charts[NB_CHARTS] = { 0 };
    mr_chart_init(&charts[0], chart0_bounds, mr_chart_euclidean_metric);

    mr_manifold *manifold = mr_manifold_create(3, NB_CHARTS);
    for (size_t i = 0; i < NB_CHARTS; ++i) {
        mr_manifold_set_chart(manifold, i, &charts[i]);
    }

    return manifold;
}

bool adaptive_refine(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    float norm = l2(node->x, node->y, node->z);

    return norm < 2.25f && norm > 1.75f;
}

bool is_node_in_bounds(mr_ocforest *forest, mr_chart *chart, mr_octree_node *node) {
    mr_float center[3] = { node->x, node->y, node->z };

    if (!chart->bounds(forest->manifold, chart, center)) {
        return false;
    }

    for (size_t i = 0; i < 8; ++i) {
        mr_float hdim = node->dim / 2.0f;
        mr_float x = center[0] + hdim * (i & 1 ? 1.0f : -1.0f);
        mr_float y = center[1] + hdim * ((i >> 1) & 1 ? 1.0f : -1.0f);
        mr_float z = center[2] + hdim * ((i >> 2) & 1 ? 1.0f : -1.0f);

        mr_float p[3] = { x, y, z };
        if (!chart->bounds(forest->manifold, chart, p)) {
            return false;
        }
    }

    return true;
}

bool chart0_activate(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_chart *chart = mr_manifold_get_chart(forest->manifold, node->chart_idx);

    return is_node_in_bounds(forest, chart, node);
}

mr_ocforest *setup_ocforest(mr_manifold *manifold) {
    mr_octree_root_desc descs[] = {
        (mr_octree_root_desc) {
            .flags = 0,
            .chart_idx = 0,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 8.0f,
        },
    };
    mr_ocforest *forest = mr_ocforest_create(manifold, descs, 1, NULL, 0);

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);


    for (size_t i = 0; i < 3; ++i) {
        mr_octree_refine_all(forest, 0);
    }
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_make_octree_cond_cb(adaptive_refine, NULL), true);
    mr_octree_balance(forest, 0);


    clock_gettime(CLOCK_MONOTONIC, &end);
    long long elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                           (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Elapsed time: %lld microseconds\n", elapsed_us);

    mr_octree_activate(forest, 0, mr_make_octree_cond_cb(chart0_activate, NULL));

    printf("Ocforest size: %lu\n", mr_ocforest_size(forest));

    return forest;
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

    mr_storage_buffer *octree_buffer = mr_storage_buffer_create(0);
    mr_storage_buffer_alloc(octree_buffer, mr_ocforest_size(forest) * sizeof(mr_octree_node),
        MR_STATIC_DRAW, mr_mem_pool_array_ptr(forest->nodes, MR_OCTREE_NODE_FIELD));

    update_userdata update_data = { window, camera, camera_buffer };
    
    mr_shader_source compute_shader_source = mr_shader_source_read("../shaders/octree_slice.comp", MR_COMPUTE);
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