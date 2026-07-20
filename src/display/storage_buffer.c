#include <stdlib.h>
#include "glad/gl.h"

#include "maniray/display/storage_buffer.h"
#include "maniray/utils/xmalloc.h"

struct mr_storage_buffer {
    size_t cap;
    GLuint buffer;
};

mr_storage_buffer *mr_storage_buffer_create(mr_uint binding) {
    mr_storage_buffer *storage_buffer = xmalloc(sizeof(mr_storage_buffer));

    storage_buffer->cap = 0;
    glGenBuffers(1, &storage_buffer->buffer);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, binding, storage_buffer->buffer);

    return storage_buffer;
}

void mr_storage_buffer_destroy(mr_storage_buffer *storage_buffer) {
    if (!storage_buffer) {
        return;
    }

    glDeleteBuffers(1, &storage_buffer->buffer);
    free(storage_buffer);
}

size_t mr_storage_buffer_get_capacity(mr_storage_buffer *storage_buffer) {
    return storage_buffer ? storage_buffer->cap : 0;
}

void mr_storage_buffer_alloc(mr_storage_buffer *storage_buffer, mr_isize size, mr_buffer_type type, const void *data) {
    if (!storage_buffer) {
        return;
    }

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, storage_buffer->buffer);
    glBufferData(GL_SHADER_STORAGE_BUFFER, size, data, BUFFER_TYPE_CONVERT[type]);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    storage_buffer->cap = (size_t)size;
}

void mr_storage_buffer_fill(mr_storage_buffer *storage_buffer, mr_isize offset, mr_isize size, const void *data) {
    if (!storage_buffer) {
        return;
    }

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, storage_buffer->buffer);
    glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset, size, data);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

void *mr_storage_buffer_map(mr_storage_buffer *storage_buffer, mr_isize offset, mr_isize size, mr_uint flags) {
    if (!storage_buffer) {
        return NULL;
    }

    GLbitfield access = 0;

    access |= flags & MR_STORAGE_BUFFER_FLAG_READ ? GL_MAP_READ_BIT : 0;
    access |= flags & MR_STORAGE_BUFFER_FLAG_WRITE ? GL_MAP_WRITE_BIT : 0;

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, storage_buffer->buffer);
    void *ptr = glMapBufferRange(GL_SHADER_STORAGE_BUFFER, offset, size, access);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    return ptr;
}

void mr_storage_buffer_unmap(mr_storage_buffer *storage_buffer) {
    if (!storage_buffer) {
        return;
    }

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, storage_buffer->buffer);
    glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}