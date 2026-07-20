#ifndef _MR_DISPLAY_MISC_H
#define _MR_DISPLAY_MISC_H

#include "glad/gl.h"

typedef enum mr_buffer_type {
    MR_STATIC_DRAW = 0,
    MR_DYNAMIC_DRAW,
} mr_buffer_type;

static const GLenum BUFFER_TYPE_CONVERT[] = {
    [MR_STATIC_DRAW] = GL_STATIC_DRAW,
    [MR_DYNAMIC_DRAW] = GL_DYNAMIC_DRAW,
};

#endif // _MR_DISPLAY_MISC_H