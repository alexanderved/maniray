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
#include "maniray/compute/octree.h"
#include "maniray/compute/sparse_matrix.h"
#include "maniray/compute/fvm_general.h"

static mr_float proj_radius(const mr_float *p, const mr_float *c) {
    return mr_norm2_2d(p[0] - c[0], p[2] - c[2]);
}

static mr_float dist_to_circle(const mr_float *p, const mr_float *c) {
    mr_float proj = proj_radius(p, c) - 1.0;
    mr_float n = p[1] - c[1];

    return mr_norm2_2d(proj, n);
}

static bool chart_0_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);

    return p[1] > 0.25f || p[1] < -0.25f || (mr_norm2_2d(p[0] + 2.0f, p[2]) > 1.0f && mr_norm2_2d(p[0] - 2.0f, p[2]) > 1.0f);
}

static bool transition_0_1_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return (p[1] < -0.25f && p[1] > -1.0f && mr_norm2_2d(p[0] + 2.0f, p[2]) < 1.0f)
        || (p[1] > 0.25f && p[1] < 1.0f && mr_norm2_2d(p[0] - 2.0f, p[2]) < 1.0f);
}

static int transition_0_1(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[1] > 0.0f) {
        p_out[0] -= 2.0f;
    } else {
        p_out[0] += 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_0_2_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return (p[1] < -0.25f && p[1] > -1.0f && mr_norm2_2d(p[0] - 2.0f, p[2]) < 1.0f)
        || (p[1] > 0.25f && p[1] < 1.0f && mr_norm2_2d(p[0] + 2.0f, p[2]) < 1.0f);
}

static int transition_0_2(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[1] > 0.0f) {
        p_out[0] += 2.0f;
    } else {
        p_out[0] -= 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_0_3_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    mr_float mouth_circle1[3] = { -2.0f, 0.0f, 0.0f };
    mr_float mouth_circle2[3] = { 2.0f, 0.0f, 0.0f };

    return dist_to_circle(p, mouth_circle1) < 0.75 || dist_to_circle(p, mouth_circle2) < 0.75;
}

static int transition_0_3(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float mouth_circle1[3] = { -2.0f, 0.0f, 0.0f };
    mr_float mouth_circle2[3] = { 2.0f, 0.0f, 0.0f };

    mr_float r1 = dist_to_circle(p_in, mouth_circle1);
    mr_float r2 = dist_to_circle(p_in, mouth_circle2);

    mr_float r = 0.0f;
    const mr_float *mc = NULL;
    if (r1 < 0.75f) {
        r = r1;
        mc = mouth_circle1;
    } else if (r2 < 0.75f) {
        r = r2;
        mc = mouth_circle2;
    } else {
        return MR_FAILURE;
    }

    mr_float theta = mr_atan2p(p_in[1] - mc[1], proj_radius(p_in, mc) - 1.0f);
    mr_float phi = mr_atan2p(p_in[2] - mc[2], p_in[0] - mc[0]);

    if ((r1 < 0.75f && theta > MR_PI) || r2 < 0.75f) {
        theta += 2.0f * MR_PI;
    }

    p_out[0] = r;
    p_out[1] = theta / (8.0f * MR_PI) - 0.25f;
    p_out[2] = phi / (4.0f * MR_PI) - 0.25f;

    return MR_SUCCESS;
}

static bool chart_1_2_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);

    return p[1] < 1.0f && p[1] > -1.0f && mr_norm2_2d(p[0], p[2]) < 1.0f;
}

static bool transition_1_0_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return p[1] < -0.25f || p[1] > 0.25f;
}

static int transition_1_0(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[1] > 0.0f) {
        p_out[0] += 2.0f;
    } else {
        p_out[0] -= 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_1_3_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return dist_to_circle(p, (mr_float[3]) { 0.0f, 0.0f, 0.0f }) < 0.75;
}

static int transition_1_3(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float mouth_circle[3] = { 0.0f, 0.0f, 0.0f };

    mr_float r = dist_to_circle(p_in, mouth_circle);
    mr_float theta = mr_atan2p(p_in[1] - mouth_circle[1], proj_radius(p_in, mouth_circle) - 1.0f) + 2.0f * MR_PI;
    mr_float phi = mr_atan2p(p_in[2] - mouth_circle[2], p_in[0] - mouth_circle[0]);

    p_out[0] = r;
    p_out[1] = theta / (8.0f * MR_PI) - 0.25f;
    p_out[2] = phi / (4.0f * MR_PI) - 0.25f;

    return MR_SUCCESS;
}

static bool transition_2_0_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return p[1] < -0.25f || p[1] > 0.25f;
}

static int transition_2_0(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[1] > 0.0f) {
        p_out[0] -= 2.0f;
    } else {
        p_out[0] += 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_2_3_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return dist_to_circle(p, (mr_float[3]) { 0.0f, 0.0f, 0.0f }) < 0.75;
}

static int transition_2_3(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float mouth_circle[3] = { 0.0f, 0.0f, 0.0f };

    mr_float r = dist_to_circle(p_in, mouth_circle);
    mr_float theta = mr_atan2p(p_in[1] - mouth_circle[1], proj_radius(p_in, mouth_circle) - 1.0f);
    mr_float phi = mr_atan2p(p_in[2] - mouth_circle[2], p_in[0] - mouth_circle[0]);

    p_out[0] = r;
    p_out[1] = theta / (8.0f * MR_PI) - 0.25f;
    p_out[2] = phi / (4.0f * MR_PI) - 0.25f;

    return MR_SUCCESS;
}

static bool chart_3_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);

    return p[0] > 0.0 && p[0] < 0.75;
}

static void chart_3_period(const mr_chart *chart, mr_float *wp, const mr_float *p) {
    MR_UNUSED(chart);

    wp[0] = p[0];
    wp[1] = mr_wrap(p[1], -0.25, 0.25);
    wp[2] = mr_wrap(p[2], -0.25, 0.25);
}

static bool transition_3_0_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    mr_float theta = (p[1] + 0.25f) * 8.0f * MR_PI;
    return theta != MR_PI && theta != 3.0f * MR_PI;
}

static int transition_3_0(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float r = p_in[0];
    mr_float theta = (p_in[1] + 0.25f) * 8.0f * MR_PI;
    mr_float phi = (p_in[2] + 0.25f) * 4.0f * MR_PI;

    mr_float x_offset = 0.0f;
    if ((theta >= 0.0f && theta < 1.0f * MR_PI) || (theta > 3.0f * MR_PI && theta <= 4.0f * MR_PI)) {
        x_offset = -2.0f;
    } else if (theta > 1.0f * MR_PI && theta < 3.0f * MR_PI) {
        x_offset = 2.0f;
    } else {
        return MR_FAILURE;
    }

    p_out[0] = (1.0f + r * cosf(theta)) * cosf(phi) + x_offset;
    p_out[1] = r * sinf(theta);
    p_out[2] = (1.0f + r * cosf(theta)) * sinf(phi);

    return MR_SUCCESS;
}

static bool transition_3_1_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    mr_float theta = (p[1] + 0.25f) * 8.0f * MR_PI;
    return theta > 2.5f * MR_PI && theta < 3.5f * MR_PI;
}

static int transition_3_1(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float r = p_in[0];
    mr_float theta = (p_in[1] + 0.25f) * 8.0f * MR_PI;
    mr_float phi = (p_in[2] + 0.25f) * 4.0f * MR_PI;

    p_out[0] = (1.0f + r * cosf(theta)) * cosf(phi);
    p_out[1] = r * sinf(theta);
    p_out[2] = (1.0f + r * cosf(theta)) * sinf(phi);

    return MR_SUCCESS;
}

static bool transition_3_2_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    mr_float theta = (p[1] + 0.25f) * 8.0f * MR_PI;
    return theta > 0.5f * MR_PI && theta < 1.5f * MR_PI;
}

static int transition_3_2(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float r = p_in[0];
    mr_float theta = (p_in[1] + 0.25f) * 8.0f * MR_PI;
    mr_float phi = (p_in[2] + 0.25f) * 4.0f * MR_PI;

    p_out[0] = (1.0f + r * cosf(theta)) * cosf(phi);
    p_out[1] = r * sinf(theta);
    p_out[2] = (1.0f + r * cosf(theta)) * sinf(phi);

    return MR_SUCCESS;
}

mr_manifold *setup_manifold() {
#define NB_CHARTS 4
    mr_chart_desc charts[NB_CHARTS] = {
        { .bounds = chart_0_bounds },
        { .bounds = chart_1_2_bounds },
        { .bounds = chart_1_2_bounds },
        { .bounds = chart_3_bounds, .period = chart_3_period }, // TODO: Add metric
    };

    mr_transition_desc transitions[NB_CHARTS * NB_CHARTS] = {
        mr_transition_desc_create_self(),
        { transition_0_1_domain, transition_0_1 },
        { transition_0_2_domain, transition_0_2 },
        { transition_0_3_domain, transition_0_3 },

        { transition_1_0_domain, transition_1_0 },
        mr_transition_desc_create_self(),
        mr_transition_desc_create_empty(),
        { transition_1_3_domain, transition_1_3 },

        { transition_2_0_domain, transition_2_0 },
        mr_transition_desc_create_empty(),
        mr_transition_desc_create_self(),
        { transition_2_3_domain, transition_2_3 },

        { transition_3_0_domain, transition_3_0 },
        { transition_3_1_domain, transition_3_1 },
        { transition_3_2_domain, transition_3_2 },
        mr_transition_desc_create_self(),
    };

    return mr_manifold_create(3, NB_CHARTS, charts, transitions);
}

static bool adaptive_refine(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_octree_node *node = mr_ocforest_get_node(forest, idx);
    mr_float p[3] = { node->x, node->y, node->z };

    return p[1] < 1.25f && p[1] > -1.25f && (mr_norm2_2d(p[0] + 2.0f, p[2]) < 1.8f || mr_norm2_2d(p[0] - 2.0f, p[2]) < 1.8f);
}

static int setup_boundary(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    for (mr_octree_direction dir = MR_OCTREE_DIRECTION_MI_X; dir <= MR_OCTREE_DIRECTION_PL_Z; ++dir) {
        mr_int nidx = mr_octree_find_face_neighbor(forest, idx, dir);
        mr_octree_node *n = mr_ocforest_get_node(forest, nidx);
        if (!n /* || !(n->flags & MR_OCTREE_NODE_FLAG_ACTIVE) */) {
            discr_data->type = MR_CELL_TYPE_BOUNDARY;

            break;
        }
    }

    return MR_SUCCESS;
}

static int setup_boundary_chart3(mr_ocforest *forest, mr_int idx, void *userdata) {
    MR_UNUSED(userdata);

    mr_discretization_data *discr_data = mr_ocforest_get_extra(forest, idx, MR_DISCR_DATA_EXTRA_FIELD);
    if (discr_data->type != MR_CELL_TYPE_NONE) {
        return MR_SUCCESS;
    }

    mr_octree_direction dir = MR_OCTREE_DIRECTION_MI_X;

    mr_int nidx = mr_octree_find_face_neighbor(forest, idx, dir);
    mr_octree_node *n = mr_ocforest_get_node(forest, nidx);

    if (!n) {
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


    for (size_t i = 0; i < 4; ++i) {
        mr_octree_refine_all(forest, 0);
        mr_octree_refine_all(forest, 1);
        mr_octree_refine_all(forest, 2);
        mr_octree_refine_all(forest, 3);
    }
    mr_octree_refine(forest, 0, mr_octree_cond_cb_null(), mr_make_octree_cond_cb(adaptive_refine, NULL), true);
    mr_octree_balance(forest, 0);


    clock_gettime(CLOCK_MONOTONIC, &end);
    long long elapsed_us = (end.tv_sec - start.tv_sec) * 1000000LL + 
                           (end.tv_nsec - start.tv_nsec) / 1000;

    printf("Refine + Balance: %.2f ms\n", (double)elapsed_us / 1000.0);


    clock_gettime(CLOCK_MONOTONIC, &start);

    mr_fvm_fit_grids_to_charts(forest);
    mr_octree_leaves_apply(forest, 0, mr_octree_cond_cb_null(), mr_make_octree_apply_cb(setup_boundary, NULL), false);
    mr_octree_leaves_apply(forest, 3, mr_octree_cond_cb_null(), mr_make_octree_apply_cb(setup_boundary_chart3, NULL), false);
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