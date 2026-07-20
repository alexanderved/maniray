#include <stdbool.h>
#include <GLFW/glfw3.h>

#include "maniray/display/camera.h"
#include "maniray/display/window_internal.h"
#include "maniray/utils/misc.h"
#include "maniray/utils/xmalloc.h"

#define SENSETIVITY 0.75f
#define HSPEED 6.5f
#define VSPEED 4.5f

struct mr_camera {
    mr_window *window;

    vec3 position;

    vec3 world_up;
    vec3 forward;
    vec3 right;
    vec3 up;
    
    float yaw;
    float pitch;

    float fov;
    float focal_length;

    bool first_movement;
    float last_x;
    float last_y;
};

static void update_camera_vectors(mr_camera *camera) {
    camera->forward[0] = sin(glm_rad(camera->yaw)) * cos(glm_rad(camera->pitch));
    camera->forward[1] = -sin(glm_rad(camera->pitch));
    camera->forward[2] = cos(glm_rad(camera->yaw)) * cos(glm_rad(camera->pitch));

    glm_vec3_normalize(camera->forward);
    glm_vec3_crossn(camera->world_up, camera->forward, camera->right);
    glm_vec3_crossn(camera->forward, camera->right, camera->up);
}

static void mouse_movement_callback(mr_window *window, double xposd, double yposd, void *userdata) {
    MR_UNUSED(window);

    float xpos = (float)xposd;
    float ypos = (float)yposd;
    mr_camera *camera = userdata;

    if (camera->first_movement) {
        camera->last_x = xpos;
        camera->last_y = ypos;
        camera->first_movement = false;
    }

    float xoffset = xpos - camera->last_x;
    float yoffset = ypos - camera->last_y; 
    camera->last_x = xpos;
    camera->last_y = ypos;

    xoffset *= SENSETIVITY;
    yoffset *= SENSETIVITY;

    camera->yaw += xoffset;
    camera->pitch += yoffset;

    // printf("OFFSETS: %f, %f\n", -xoffset, -yoffset);

    if (camera->pitch > 89.0f) {
        camera->pitch = 89.0f;
    }

    if (camera->pitch < -89.0f) {
        camera->pitch = -89.0f;
    }

    //printf("%f, %f\n", camera->yaw, camera->pitch);

    update_camera_vectors(camera);
}

static void key_press_callback(mr_window *window, int key, int action, void *userdata) {
    mr_camera *camera = userdata;
    float delta_time = mr_window_get_delta_time(window);

    vec3 front = { camera->forward[0], 0.0f, camera->forward[2] };
    glm_vec3_normalize(front);

    if (key == GLFW_KEY_W && action == GLFW_PRESS) {
        vec3 delta = { 0 };
        glm_vec3_scale(front, HSPEED * delta_time, delta);
        glm_vec3_add(
            camera->position,
            delta,
            camera->position);
    }

    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        vec3 delta = { 0 };
        glm_vec3_scale(front, HSPEED * delta_time, delta);
        glm_vec3_sub(
            camera->position,
            delta,
            camera->position);
    }
    
    if (key == GLFW_KEY_A && action == GLFW_PRESS) {
        vec3 delta = { 0 };
        glm_vec3_scale(camera->right, HSPEED * delta_time, delta);
        glm_vec3_sub(
            camera->position,
            delta,
            camera->position);
    }

    if (key == GLFW_KEY_D && action == GLFW_PRESS) {
        vec3 delta = { 0 };
        glm_vec3_scale(camera->right, HSPEED * delta_time, delta);
        glm_vec3_add(
            camera->position,
            delta,
            camera->position);
    }

    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
        vec3 delta = { 0 };
        glm_vec3_scale(camera->world_up, VSPEED * delta_time, delta);
        glm_vec3_add(
            camera->position,
            delta,
            camera->position);
    }

    if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_PRESS) {
        vec3 delta = { 0 };
        glm_vec3_scale(camera->world_up, VSPEED * delta_time, delta);
        glm_vec3_sub(
            camera->position,
            delta,
            camera->position);
    }
}

mr_camera *mr_camera_create(mr_window *window, vec3 position, float fov, float focal_length) {
    mr_camera *camera = xmalloc(sizeof(mr_camera));

    camera->window = window;
    glm_vec3_copy(position, camera->position);

    camera->yaw = 0.0f;
    camera->pitch = 0.0f;

    glm_vec3_copy((vec3){ 0.0, 1.0, 0.0 }, camera->world_up);
    update_camera_vectors(camera);

    camera->fov = fov;
    camera->focal_length = focal_length;

    mr_window_resolution res = mr_window_get_resolution(window);
    camera->first_movement = true;
    camera->last_x = res.width / 2.0;
    camera->last_y = res.height / 2.0;

    // mr_window_set_key_callback(window, key_press_callback, camera);
    mr_window_set_mouse_callback(window, mouse_movement_callback, camera);

    return camera;
}

void mr_camera_destroy(mr_camera *camera) {
    if (!camera) {
        return;
    }

    free(camera);
}

mr_camera_info mr_camera_get_info(mr_camera *camera) {
    mat3 yaw_rot = {
        { cos(glm_rad(camera->yaw)), 0, -sin(glm_rad(camera->yaw)) },
        { 0, 1, 0 },
        { sin(glm_rad(camera->yaw)), 0, cos(glm_rad(camera->yaw)) }
    };
    mat3 pitch_rot = {
        { 1, 0, 0 },
        { 0, cos(glm_rad(camera->pitch)), sin(glm_rad(camera->pitch)) },
        { 0, -sin(glm_rad(camera->pitch)), cos(glm_rad(camera->pitch)) }
    };

    mat3 rotation = { 0 };
    glm_mat3_mul(yaw_rot, pitch_rot, rotation);

    mat3x4 rotation_aligned = { 0 };
    glm_vec4(rotation[0], 0.0f, rotation_aligned[0]);
    glm_vec4(rotation[1], 0.0f, rotation_aligned[1]);
    glm_vec4(rotation[2], 0.0f, rotation_aligned[2]);

    mr_camera_info info;
    glm_vec3_copy(camera->position, info.position);
    glm_mat3x4_copy(rotation_aligned, info.rotation);
    info.fov = camera->fov;
    info.focal_length = camera->focal_length;

    return info;
}

void mr_camera_proccess_movement(mr_camera *camera, mr_window *window) {
    if (glfwGetKey(window->glfw_window, GLFW_KEY_W) == GLFW_PRESS)
        key_press_callback(window, GLFW_KEY_W, GLFW_PRESS, camera);

    if (glfwGetKey(window->glfw_window, GLFW_KEY_S) == GLFW_PRESS)
        key_press_callback(window, GLFW_KEY_S, GLFW_PRESS, camera);

    if (glfwGetKey(window->glfw_window, GLFW_KEY_A) == GLFW_PRESS)
        key_press_callback(window, GLFW_KEY_A, GLFW_PRESS, camera);

    if (glfwGetKey(window->glfw_window, GLFW_KEY_D) == GLFW_PRESS)
        key_press_callback(window, GLFW_KEY_D, GLFW_PRESS, camera);

    if (glfwGetKey(window->glfw_window, GLFW_KEY_SPACE) == GLFW_PRESS)
        key_press_callback(window, GLFW_KEY_SPACE, GLFW_PRESS, camera);

    if (glfwGetKey(window->glfw_window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        key_press_callback(window, GLFW_KEY_LEFT_SHIFT, GLFW_PRESS, camera);
}