#ifndef _MR_CAMERA_H
#define _MR_CAMERA_H

#include "cglm/cglm.h"

#include "window.h"

typedef struct mr_camera mr_camera;

typedef struct mr_camera_info {
    vec3 position;
    float _pad;
    mat3x4 rotation;
    float fov;
    float focal_length;
} mr_camera_info;

mr_camera *mr_camera_create(mr_window *window, vec3 position, float fov, float focal_length);
void mr_camera_destroy(mr_camera *camera);

mr_camera_info mr_camera_get_info(mr_camera *camera);
void mr_camera_proccess_movement(mr_camera *camera, mr_window *window);

#endif // _MR_CAMERA_H