#ifndef _MR_WINDOW_INTERNAL_H
#define _MR_WINDOW_INTERNAL_H

#include <GLFW/glfw3.h>

#include "engine.h"

typedef struct mr_window_userdata {
    mr_window *window;
    mr_engine *engine;

    mr_key_callback key_callback_func;
    void *key_callback_userdata;

    mr_mouse_callback mouse_callback_func;
    void *mouse_callback_userdata;
} mr_window_userdata;

struct mr_window {
    GLFWwindow *glfw_window;

    mr_window_userdata userdata;

    float last_frame;
    float delta_time;
};

#endif // _MR_WINDOW_INTERNAL_H