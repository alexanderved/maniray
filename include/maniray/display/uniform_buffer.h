#ifndef _MR_UNIFORM_BUFFER_H
#define _MR_UNIFORM_BUFFER_H

#include "maniray/utils/types.h"
#include "maniray/display/misc.h"

typedef struct mr_uniform_buffer mr_uniform_buffer;

mr_uniform_buffer *mr_uniform_buffer_create(mr_uint binding, mr_isize size, mr_buffer_type type);
void mr_uniform_buffer_destroy(mr_uniform_buffer *uniform_buffer);

void mr_uniform_buffer_fill(mr_uniform_buffer *uniform_buffer, mr_isize offset, mr_isize size, const void *data);

#endif // _MR_UNIFORM_BUFFER_H