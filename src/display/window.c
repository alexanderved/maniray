#include <stdlib.h>

#include "maniray/display/window_internal.h"
#include "maniray/utils/misc.h"
#include "maniray/utils/xmalloc.h"

void mr_window_init() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
}

void mr_window_terminate() {
    glfwTerminate();
}

static void default_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    MR_UNUSED(scancode);
    MR_UNUSED(mods);

    mr_window_userdata *userdata = glfwGetWindowUserPointer(window);
    if (userdata->key_callback_func) {
        userdata->key_callback_func(userdata->window, key, action, userdata->key_callback_userdata);
    }

    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
}

static void default_mouse_callback(GLFWwindow *window, double xpos, double ypos) {
    mr_window_userdata *userdata = glfwGetWindowUserPointer(window);
    if (userdata->mouse_callback_func) {
        userdata->mouse_callback_func(userdata->window, xpos, ypos, userdata->mouse_callback_userdata);
    }
}

mr_window *mr_window_create(int width, int height, const char *title) {
    mr_window *window = xmalloc(sizeof(mr_window));

    window->glfw_window = glfwCreateWindow(width, height, title, NULL, NULL);
    if (!window->glfw_window) {
        free(window);

        return NULL;
    }

    window->userdata = (mr_window_userdata){ 0 };
    window->userdata.window = window;
    glfwSetWindowUserPointer(window->glfw_window, &window->userdata);

    window->last_frame = 0.0f;
    window->delta_time = 0.0f;

    glfwSetKeyCallback(window->glfw_window, default_key_callback);
    glfwSetCursorPosCallback(window->glfw_window, default_mouse_callback);
    
    return window;
}

void mr_window_destroy(mr_window *window) {
    if (!window) {
        return;
    }

    glfwSetWindowUserPointer(window->glfw_window, NULL);
    if (!window->glfw_window) {
        glfwDestroyWindow(window->glfw_window);
    }
    free(window);
}

mr_window_resolution mr_window_get_resolution(mr_window *window) {
    int width, height;
    glfwGetWindowSize(window->glfw_window, &width, &height);

    return (mr_window_resolution){ width, height };
}

float mr_window_get_delta_time(mr_window *window) {
    return window->delta_time;
}

void mr_window_set_cursor_hidden(mr_window *window, bool hidden) {
    glfwSetInputMode(window->glfw_window, GLFW_CURSOR, hidden ? GLFW_CURSOR_DISABLED : GLFW_CURSOR_NORMAL);
}

void mr_window_set_key_callback(mr_window *window, mr_key_callback func, void *userdata) {
    window->userdata.key_callback_func = func;
    window->userdata.key_callback_userdata = userdata;
}

void mr_window_set_mouse_callback(mr_window *window, mr_mouse_callback func, void *userdata) {
    window->userdata.mouse_callback_func = func;
    window->userdata.mouse_callback_userdata = userdata;
}

void mr_window_poll(mr_window *window) {
    float current_frame = (float)glfwGetTime();
    window->delta_time = current_frame - window->last_frame;
    window->last_frame = current_frame;

    glfwSwapBuffers(window->glfw_window);
    glfwPollEvents();
}