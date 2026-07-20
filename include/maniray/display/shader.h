#ifndef _MR_SHADER_H
#define _MR_SHADER_H

#include <stdbool.h>

typedef enum mr_shader_type {
    MR_VERTEX = 0,
    MR_FRAGMENT,
    MR_COMPUTE,
} mr_shader_type;

typedef struct mr_shader_source {
    char *source;
    mr_shader_type type;
} mr_shader_source;

mr_shader_source mr_shader_source_create(const char *source, mr_shader_type type);
void mr_shader_source_destroy(mr_shader_source source);

mr_shader_source mr_shader_source_read(const char *filepath, mr_shader_type type);
bool mr_shader_source_is_valid(mr_shader_source source);

#endif // _MR_SHADER_H