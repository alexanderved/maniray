#include <stdlib.h>
#include <stdio.h>
#include "glad/gl.h"
#include <GLFW/glfw3.h>

#include "maniray/display/engine.h"
#include "maniray/display/window_internal.h"
#include "maniray/display/program.h"
#include "maniray/display/uniform_buffer.h"
#include "maniray/utils/misc.h"
#include "maniray/utils/xmalloc.h"

typedef struct engine_info {
    float time;
    int width;
    int height;
} engine_info;

struct mr_engine {
    mr_window *window;

    mr_program *display_program;
    GLuint vbo;
    GLuint vao;

    GLuint texture;

    mr_uniform_buffer *info_buffer;
};

int mr_engine_init(mr_window *window) {
    glfwMakeContextCurrent(window->glfw_window);
    glfwSwapInterval(0);

    return gladLoadGL((GLADloadfunc)glfwGetProcAddress);
}

static void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    mr_window_userdata *userdata = glfwGetWindowUserPointer(window);
    mr_engine *engine = userdata->engine;

    glViewport(0, 0, width, height);

    mr_window_resolution resolution = { width, height };
    mr_uniform_buffer_fill(engine->info_buffer, offsetof(engine_info, width), sizeof(resolution), &resolution);
}

#define NB_VERTICES 4
#define NB_COORDS 2

static void fill_vbo(GLuint vbo) {
    const float vertices[NB_VERTICES * NB_COORDS] = {
         1.0f,  1.0f,
         1.0f, -1.0f,
        -1.0f,  1.0f,
        -1.0f, -1.0f,
    };

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
}

static void set_attribute_pointer() {
    glVertexAttribPointer(0, NB_COORDS, GL_FLOAT, GL_FALSE, NB_COORDS * sizeof(float), NULL);
    glEnableVertexAttribArray(0);
}

static void setup_buffers(mr_engine *engine) {
    glGenVertexArrays(1, &engine->vao);
    glGenBuffers(1, &engine->vbo);

    glBindVertexArray(engine->vao);
    fill_vbo(engine->vbo);
    set_attribute_pointer();

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

#define TEXTURE_SIZE 1024

static void create_texture(mr_engine *engine) {
    glGenTextures(1, &engine->texture);
    glBindTexture(GL_TEXTURE_2D, engine->texture);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glTexImage2D(
        GL_TEXTURE_2D,
        0,
        GL_RGBA32F,
        TEXTURE_SIZE,
        TEXTURE_SIZE,
        0,
        GL_RGBA,
        GL_FLOAT,
        NULL
    );

    glBindTexture(GL_TEXTURE_2D, 0);
}

mr_engine *mr_engine_create(mr_window *window) {
    mr_engine *engine = xmalloc(sizeof(mr_engine));
    engine->window = window;

    mr_shader_source display_shaders[] = {
        mr_shader_source_read("../shaders/display.vert", MR_VERTEX),
        mr_shader_source_read("../shaders/display.frag", MR_FRAGMENT),
    };
    engine->display_program = mr_program_create(display_shaders, 2);

    mr_shader_source_destroy(display_shaders[0]);
    mr_shader_source_destroy(display_shaders[1]);

    if (!engine->display_program) {
        free(engine);

        return NULL;
    }

    setup_buffers(engine);
    create_texture(engine);

    engine->info_buffer = mr_uniform_buffer_create(
        MR_INFO_BINDING, sizeof(engine_info), MR_DYNAMIC_DRAW);
    if (!engine->info_buffer) {
        mr_program_destroy(engine->display_program);
        free(engine);

        return NULL;
    }

    window->userdata.engine = engine;
    glfwSetFramebufferSizeCallback(window->glfw_window, framebuffer_size_callback);

    return engine;
}

void mr_engine_destroy(mr_engine *engine) {
    if (!engine) {
        return;
    }

    mr_uniform_buffer_destroy(engine->info_buffer);
    glDeleteTextures(1, &engine->texture);
    glDeleteVertexArrays(1, &engine->vao);
    glDeleteBuffers(1, &engine->vbo);
    mr_program_destroy(engine->display_program);

    free(engine);
}

void mr_engine_run(mr_engine *engine, mr_program *compute_program, mr_engine_update update, void *userdata) {
    if (!engine || !engine->window || !engine->window->glfw_window) {
        return;
    }

    mr_window_resolution resolution = mr_window_get_resolution(engine->window);
    engine_info info = {
        .time = 0.0f,
        .width = resolution.width,
        .height = resolution.height,
    };
    mr_uniform_buffer_fill(engine->info_buffer, 0, sizeof(info), &info);

    glViewport(0, 0, resolution.width, resolution.height);

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    glBindTexture(GL_TEXTURE_2D, engine->texture);
    glBindImageTexture(MR_TEXTURE_BINDING, engine->texture, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);

    while (!glfwWindowShouldClose(engine->window->glfw_window)) {
        float time = (float)glfwGetTime();
        mr_uniform_buffer_fill(engine->info_buffer, offsetof(engine_info, time), sizeof(float), &time);

        mr_program_use(compute_program);
        glDispatchCompute(TEXTURE_SIZE / 8, TEXTURE_SIZE / 8, 1);

        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

        glClearColor(0.9f, 0.9f, 0.9f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        mr_program_use(engine->display_program);
        glBindVertexArray(engine->vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, NB_VERTICES);

        update(engine, userdata);
        mr_window_poll(engine->window);
    }
}