#include <stdio.h>
#include <time.h>

#include "maniray/display/engine.h"
#include "maniray/display/camera.h"
#include "maniray/display/uniform_buffer.h"
#include "maniray/display/storage_buffer.h"
#include "maniray/utils/misc.h"
#include "maniray/compute/octree.h"

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

float l2(float x, float y, float z) {
    return sqrtf(x * x + y * y + z * z);
}

float linf(float x, float y, float z) {
    return MAX(MAX(fabsf(x), fabsf(y)), fabsf(z));
}

float cube_sdf(mr_ocforest *forest, mr_int idx) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    float x = node->x;
    float y = node->y;
    float z = node->z;

    float dim = 0.5;
    float ax = fabsf(x) - dim;
    float ay = fabsf(y) - dim;
    float az = fabsf(z) - dim;

    float mx = MAX(ax, 0.0);
    float my = MAX(ay, 0.0);
    float mz = MAX(az, 0.0);
    float l = l2(mx, my, mz);

    return l + MIN(MAX(ax, MAX(ay, az)), 0.0);
}

int write_dist(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    float x = node->x;
    float y = node->y;
    float z = node->z;

    float ax = fabsf(x);
    float ay = fabsf(y);
    float az = fabsf(z);

    float norm = l2(x, y, z);
    float sdf = norm - 1.0f;
    // float sdf = linf(x, y, z) - 0.5f;
    // float sdf = cube_sdf(forest, idx);

    node->value = sdf;
    node->gradient[0] = x / norm;
    node->gradient[1] = y / norm;
    node->gradient[2] = z / norm;

    /* float normi = linf(x, y, z);
    if (normi == ax) {
        node->gradient[0] = copysignf(1.0f, x);
    }

    if (normi == ay) {
        node->gradient[1] = copysignf(1.0f, y);
    }

    if (normi == az) {
        node->gradient[2] = copysignf(1.0f, z);
    } */

    /* node->gradient[0] /= l2(node->gradient[0], node->gradient[1], node->gradient[2]);
    node->gradient[1] /= l2(node->gradient[0], node->gradient[1], node->gradient[2]);
    node->gradient[2] /= l2(node->gradient[0], node->gradient[1], node->gradient[2]); */

    return MR_SUCCESS;
}

mr_int found_node = 0;

bool adaptive_refine(mr_ocforest *forest, mr_int idx, void *userdata) {
    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    float norm = l2(node->x, node->y, node->z);

    return norm < 2.25f && norm > 1.75f;
}

#define DIM 8.0f

bool point_refine(mr_ocforest *forest, mr_int idx, void *userdata) {
    float mul = DIM;
    float x = 0.32f * mul;
    float y = 0.07f * mul;
    float z = -0.24f * mul;

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);

    if ((node->flags & MR_OCTREE_NODE_FLAG_LEAF) && (linf(x - node->x, y - node->y, z - node->z) <= node->dim / 2.0f)) {
        found_node = idx;
    }

    return linf(x - node->x, y - node->y, z - node->z) <= node->dim / 2.0f;
}

mr_ocforest *setup_ocforest() {
    mr_octree_root_desc descs[] = {
        (mr_octree_root_desc) {
            .flags = 0,
            .chart_idx = 0,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = DIM,
        },
    };
    mr_ocforest *forest = mr_ocforest_create(NULL, descs, 1);

    /* for (size_t i = 0; i < MR_OCTREE_MAX_LEVEL - 1; ++i) {
        mr_octree_refine_all(forest, 0);
    }
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_make_octree_cond_cb(adaptive_refine, NULL)); */

    /* for (size_t i = 0; i < MR_OCTREE_MAX_LEVEL; ++i) {
        mr_octree_refine(forest, 0, mr_make_octree_cond_cb(point_refine, NULL), mr_octree_cond_cb_null());
    } */

    // mr_octree_refine(forest, 0, mr_make_octree_cond_cb(point_refine, NULL), mr_octree_cond_cb_null(), true);


    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);


    mr_octree_refine_all(forest, 0);
    mr_octree_refine_all(forest, 0);
    mr_octree_refine_all(forest, 0);
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_make_octree_cond_cb(adaptive_refine, NULL), true);

    mr_octree_balance(forest, 0);
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_make_octree_cond_cb(adaptive_refine, NULL), true);
    mr_octree_balance(forest, 0);


    clock_gettime(CLOCK_MONOTONIC, &end);
    long long elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                           (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Elapsed time: %lld microseconds\n", elapsed_us);



    mr_octree_leaves_apply(forest, 0, mr_octree_cond_cb_null(), mr_make_octree_apply_cb(write_dist, NULL), false);

    mr_int node_idx = mr_octree_find_face_neighbor(forest, found_node, MR_OCTREE_DIRECTION_MI_Y);
    mr_octree_node *node = mr_ocforest_get_node(forest, node_idx);

    // printf("%u: %d\n", node->level, node_idx);

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


    mr_ocforest *forest = setup_ocforest();

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

int test_mem_pool() {
    typedef struct example {
        int a;
        float b;
    } example;

    mr_mem_pool *pool = mr_mem_pool_create((size_t[]) { sizeof(example), sizeof(char) }, 2);

    mr_mem_pool_alloc(pool);

    char c = 'a';
    mr_mem_pool_insert(pool, &(example) { 1, 2.0f }, &c);

    mr_mem_pool_alloc_many(pool, 13);
    mr_mem_pool_alloc_many(pool, 43);

    mr_mem_pool_remove(pool, 5);
    mr_mem_pool_remove(pool, 7);

    return 0;
}

int test_octree() {
    mr_octree_root_desc descs[] = {
        (mr_octree_root_desc) {
            .flags = 0,
            .chart_idx = 0,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 1.0f,
        },
    };
    mr_ocforest *forest = mr_ocforest_create(NULL, descs, 1);

    for (size_t i = 0; i < 3; ++i) {
        mr_octree_refine_all(forest, 0);
    }

    mr_ocforest_destroy(forest);

    return 0;
}

int main() {
    run_display();
    // test_mem_pool();
    // test_octree();

    return 0;
}