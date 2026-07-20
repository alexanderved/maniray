#ifndef _MR_STORAGE_BUFFER_H
#define _MR_STORAGE_BUFFER_H

#include "maniray/utils/types.h"
#include "maniray/display/misc.h"

#define MR_STORAGE_BUFFER_FLAG_READ (1 << 0)
#define MR_STORAGE_BUFFER_FLAG_WRITE (1 << 1)

typedef struct mr_storage_buffer mr_storage_buffer;

mr_storage_buffer *mr_storage_buffer_create(mr_uint binding);
void mr_storage_buffer_destroy(mr_storage_buffer *storage_buffer);

size_t mr_storage_buffer_get_capacity(mr_storage_buffer *storage_buffer);

void mr_storage_buffer_alloc(mr_storage_buffer *storage_buffer, mr_isize size, mr_buffer_type type, const void *data);
void mr_storage_buffer_fill(mr_storage_buffer *storage_buffer, mr_isize offset, mr_isize size, const void *data);

void *mr_storage_buffer_map(mr_storage_buffer *storage_buffer, mr_isize offset, mr_isize size, mr_uint flags);
void mr_storage_buffer_unmap(mr_storage_buffer *storage_buffer);

#endif // _MR_STORAGE_BUFFER_H