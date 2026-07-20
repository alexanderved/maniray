#ifndef _MR_PROGRAM_H
#define _MR_PROGRAM_H

#include "shader.h"

typedef struct mr_program mr_program;

mr_program *mr_program_create(mr_shader_source *shaders, int nb_shaders);
void mr_program_destroy(mr_program *program);

void mr_program_use(mr_program *program);

#endif // _MR_PROGRAM_H
