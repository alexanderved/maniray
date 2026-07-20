#include <stdlib.h>
#include "glad/gl.h"

#include "maniray/display/uniform_buffer.h"
#include "maniray/utils/xmalloc.h"

struct mr_uniform_buffer {
    GLuint buffer;
};

mr_uniform_buffer *mr_uniform_buffer_create(mr_uint binding, mr_isize size, mr_buffer_type type) {
    mr_uniform_buffer *uniform_buffer = xmalloc(sizeof(mr_uniform_buffer));

    glGenBuffers(1, &uniform_buffer->buffer);

    glBindBuffer(GL_UNIFORM_BUFFER, uniform_buffer->buffer);
    glBufferData(GL_UNIFORM_BUFFER, size, NULL, BUFFER_TYPE_CONVERT[type]);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

    glBindBufferBase(GL_UNIFORM_BUFFER, binding, uniform_buffer->buffer);

    return uniform_buffer;
}

void mr_uniform_buffer_destroy(mr_uniform_buffer *uniform_buffer) {
    if (!uniform_buffer) {
        return;
    }

    glDeleteBuffers(1, &uniform_buffer->buffer);
    free(uniform_buffer);
}

void mr_uniform_buffer_fill(mr_uniform_buffer *uniform_buffer, mr_isize offset, mr_isize size, const void *data) {
    glBindBuffer(GL_UNIFORM_BUFFER, uniform_buffer->buffer);
    glBufferSubData(GL_UNIFORM_BUFFER, offset, size, data);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}