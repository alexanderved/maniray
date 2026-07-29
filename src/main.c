#include <stdio.h>
#include <string.h>
#include <time.h>

#include "maniray/display/engine.h"
#include "maniray/display/camera.h"
#include "maniray/display/uniform_buffer.h"
#include "maniray/display/storage_buffer.h"
#include "maniray/utils/misc.h"
#include "maniray/compute/octree.h"
#include "maniray/compute/sparse_matrix.h"

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

static bool chart_0_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);

    return p[1] > 0.25f || p[1] < -0.25f || (l2_2d(p[0] + 2.0f, p[2]) > 1.0f && l2_2d(p[0] - 2.0f, p[2]) > 1.0f);
}

static bool chart_1_2_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);

    return p[1] < 1.0f && p[1] > -1.0f && l2_2d(p[0], p[2]) < 1.0f;
}

static bool transition_0_1_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return (p[1] < -0.25f && p[1] > -1.0f && l2_2d(p[0] + 2.0f, p[2]) < 1.0f)
        || (p[1] > 0.25f && p[1] < 1.0f && l2_2d(p[0] - 2.0f, p[2]) < 1.0f);
}

static int transition_0_1(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    if (!t->domain(t, p_in)) {
        return MR_FAILURE;
    }

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[0] > 0.0f) {
        p_out[0] -= 2.0f;
    } else {
        p_out[0] += 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_1_0_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return p[1] < -0.25f || p[1] > 0.25f;
}

static int transition_1_0(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    if (!t->domain(t, p_in)) {
        return MR_FAILURE;
    }

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[0] > 0.0f) {
        p_out[0] += 2.0f;
    } else {
        p_out[0] -= 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_0_2_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return (p[1] < -0.25f && p[1] > -1.0f && l2_2d(p[0] - 2.0f, p[2]) < 1.0f)
        || (p[1] > 0.25f && p[1] < 1.0f && l2_2d(p[0] + 2.0f, p[2]) < 1.0f);
}

static int transition_0_2(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    if (!t->domain(t, p_in)) {
        return MR_FAILURE;
    }

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[0] > 0.0f) {
        p_out[0] += 2.0f;
    } else {
        p_out[0] -= 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_2_0_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return p[1] < -0.25f || p[1] > 0.25f;
}

static int transition_2_0(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    if (!t->domain(t, p_in)) {
        return MR_FAILURE;
    }

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[0] > 0.0f) {
        p_out[0] -= 2.0f;
    } else {
        p_out[0] += 2.0f;
    }

    return MR_SUCCESS;
}

mr_manifold *setup_manifold() {
#define NB_CHARTS 3
    mr_chart_desc charts[NB_CHARTS] = {
        { .bounds = chart_0_bounds },
        { .bounds = chart_1_2_bounds },
        { .bounds = chart_1_2_bounds },
    };

    mr_transition_desc transitions[NB_CHARTS * NB_CHARTS] = {
        mr_transition_desc_create_self(), { transition_0_1_domain, transition_0_1 }, { transition_0_2_domain, transition_0_2 },
        { transition_1_0_domain, transition_1_0 }, mr_transition_desc_create_self(), mr_transition_desc_create_empty(),
        { transition_2_0_domain, transition_2_0 }, mr_transition_desc_create_empty(), mr_transition_desc_create_self(),
    };

    return mr_manifold_create(3, NB_CHARTS, charts, transitions);
}

bool adaptive_refine(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_float p[3] = { node->x, node->y, node->z };

    return p[1] < 1.25f && p[1] > -1.25f && (l2_2d(p[0] + 2.0f, p[2]) < 1.5f || l2_2d(p[0] - 2.0f, p[2]) < 1.5f);
}

bool is_node_in_bounds(mr_ocforest *forest, size_t chart_idx, mr_octree_node *node) {
    mr_float center[3] = { node->x, node->y, node->z };

    if (!mr_manifold_is_in_bounds(forest->manifold, chart_idx, center)) {
        return false;
    }

    for (size_t i = 0; i < 8; ++i) {
        mr_float hdim = node->dim / 2.0f;
        mr_float x = center[0] + hdim * (i & 1 ? 1.0f : -1.0f);
        mr_float y = center[1] + hdim * ((i >> 1) & 1 ? 1.0f : -1.0f);
        mr_float z = center[2] + hdim * ((i >> 2) & 1 ? 1.0f : -1.0f);

        mr_float p[3] = { x, y, z };
        if (!mr_manifold_is_in_bounds(forest->manifold, chart_idx, p)) {
            return false;
        }
    }

    return true;
}

bool chart_activate(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    return is_node_in_bounds(forest, node->chart_idx, node);
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
        (mr_octree_root_desc) {
            .flags = 0,
            .chart_idx = 1,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 2.0f,
        },
        (mr_octree_root_desc) {
            .flags = 0,
            .chart_idx = 2,
            .x = 0.0f,
            .y = 0.0f,
            .z = 0.0f,
            .dim = 2.0f,
        },
    };
    mr_ocforest *forest = mr_ocforest_create(manifold, descs, 3, NULL, 0);

    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);


    for (size_t i = 0; i < 4; ++i) {
        mr_octree_refine_all(forest, 0);
        mr_octree_refine_all(forest, 1);
        mr_octree_refine_all(forest, 2);
    }
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_make_octree_cond_cb(adaptive_refine, NULL), true);
    mr_octree_balance(forest, 0);


    clock_gettime(CLOCK_MONOTONIC, &end);
    long long elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                           (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Elapsed time: %lld microseconds\n", elapsed_us);

    mr_octree_activate(forest, 0, mr_make_octree_cond_cb(chart_activate, NULL));
    mr_octree_activate(forest, 1, mr_make_octree_cond_cb(chart_activate, NULL));
    mr_octree_activate(forest, 2, mr_make_octree_cond_cb(chart_activate, NULL));

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

    printf("FPS: %f\n", 1.0 / mr_window_get_delta_time(window));
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