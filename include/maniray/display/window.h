#ifndef _MR_WINDOW_H
#define _MR_WINDOW_H

#include <stdbool.h>

#define MR_WINDOW_DEFAULT_WIDTH 800
#define MR_WINDOW_DEFAULT_HEIGHT 800

typedef struct mr_window mr_window;

typedef void (*mr_key_callback)(mr_window *, int, int, void *);
typedef void (*mr_mouse_callback)(mr_window *, double, double, void *);

typedef struct mr_window_resolution {
    int width;
    int height;
} mr_window_resolution;

void mr_window_init();
void mr_window_terminate();

mr_window *mr_window_create(int width, int height, const char *title);
void mr_window_destroy(mr_window *window);

mr_window_resolution mr_window_get_resolution(mr_window *window);
float mr_window_get_delta_time(mr_window *window);

void mr_window_set_cursor_hidden(mr_window *window, bool hidden);

void mr_window_set_key_callback(mr_window *window, mr_key_callback func, void *userdata);
void mr_window_set_mouse_callback(mr_window *window, mr_mouse_callback func, void *userdata);

void mr_window_poll(mr_window *window);

#endif // _MR_WINDOW_H