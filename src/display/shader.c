#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "maniray/display/shader.h"
#include "maniray/utils/xmalloc.h"

mr_shader_source mr_shader_source_create(const char *source, mr_shader_type type) {
    char *source_copy = xmalloc(strlen(source) + 1);
    strcpy(source_copy, source);

    return (mr_shader_source){ source_copy, type };
}

void mr_shader_source_destroy(mr_shader_source source) {
    free(source.source);
}

mr_shader_source mr_shader_source_read(const char *filepath, mr_shader_type type) {
    FILE* file = fopen(filepath, "r");
    if (file == NULL) {
        return (mr_shader_source){ NULL, type };
    }

    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    rewind(file);

    char* buffer = (char*)xmalloc(length + 1);
    size_t bytes_read = fread(buffer, 1, length, file);
    if (bytes_read != (size_t)length) {
        free(buffer);
        fclose(file);

        return (mr_shader_source){ NULL, type };
    }

    buffer[length] = '\0';
    fclose(file);

    return (mr_shader_source){ buffer, type };
}

bool mr_shader_source_is_valid(mr_shader_source source) {
    return source.source != NULL;
}