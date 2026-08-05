#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "maniray/display/engine.h"
#include "maniray/display/camera.h"
#include "maniray/display/uniform_buffer.h"
#include "maniray/display/storage_buffer.h"

#include "portals_manifold.h"
#include "portals_ocforest.h"

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
    
    mr_shader_source compute_shader_source = mr_shader_source_read("../examples/portals/portals_slice.comp", MR_COMPUTE);
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