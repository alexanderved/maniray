#include <stdlib.h>
#include <stdio.h>
#include "glad/gl.h"

#include "maniray/display/program.h"
#include "maniray/utils/misc.h"
#include "maniray/utils/xmalloc.h"

struct mr_program {
    GLuint program;
};

static const GLenum SHADER_TYPE_CONVERT[] = {
    [MR_VERTEX] = GL_VERTEX_SHADER,
    [MR_FRAGMENT] = GL_FRAGMENT_SHADER,
    [MR_COMPUTE] = GL_COMPUTE_SHADER,
};

static GLuint compile_shader(const char *source, GLenum type) {
    GLuint shader = glCreateShader(type);

    glShaderSource(shader, 1, &source, NULL);
    glCompileShader(shader);

    int success;
    char info_log[512];

    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(shader, 512, NULL, info_log);
        printf("SHADER COMPILATION FAILED: %s\n", info_log);

        glDeleteShader(shader);

        return MR_FAILURE;
    }

    return shader;
}

static GLuint link_program(GLuint *shader_ids, int nb_shaders) {
    GLuint program = glCreateProgram();

    for (int i = 0; i < nb_shaders; ++i) {
        glAttachShader(program, shader_ids[i]);
    }

    glLinkProgram(program);

    int success;
    char info_log[512];

    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if(!success) {
        glGetProgramInfoLog(program, 512, NULL, info_log);
        printf("PROGRAM LINK FAILED: %s\n", info_log);

        glDeleteProgram(program);

        return MR_FAILURE;
    }

    return program;
}

mr_program *mr_program_create(mr_shader_source *shaders, int nb_shaders) {
    mr_program *program = xmalloc(sizeof(program));
    GLuint *shader_ids = xcalloc(nb_shaders, sizeof(GLuint));

    for (int i = 0; i < nb_shaders; ++i) {
        shader_ids[i] = compile_shader(shaders[i].source, SHADER_TYPE_CONVERT[shaders[i].type]);

        if (shader_ids[i] == MR_FAILURE) {
            for (int j = 0; j < i; ++j) {
                glDeleteShader(shader_ids[j]);
            }

            free(shader_ids);
            free(program);

            return NULL;
        }
    }
    program->program = link_program(shader_ids, nb_shaders);

    for (int i = 0; i < nb_shaders; ++i) {
        glDeleteShader(shader_ids[i]);
    }
    free(shader_ids);
    shader_ids = NULL;

    if (program->program == MR_FAILURE) {
        free(program);

        return NULL;
    }

    return program;
}

void mr_program_destroy(mr_program *program) {
    if (!program) {
        return;
    }

    glDeleteProgram(program->program);
    free(program);
}

void mr_program_use(mr_program *program) {
    if (!program) {
        glUseProgram(0);
    }

    glUseProgram(program->program);
}