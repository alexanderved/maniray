#ifndef _MR_DYN_ARRAY_H
#define _MR_DYN_ARRAY_H

#include <stdlib.h>

#define mr_dyn_array_create(type) mr_dyn_array_create_with_info(sizeof(type), 0, 0, NULL)
#define mr_dyn_array_create_with_capacity(type, cap) mr_dyn_array_create_with_info(sizeof(type), cap, 0, NULL)
#define mr_dyn_array_create_with_data(type, data, len) mr_dyn_array_create_with_info(sizeof(type), len, len, data)

#define mr_dyn_array_push(arr, elem) ((arr) = mr_dyn_array_push_alloc((arr), &(elem)))

void *mr_dyn_array_create_with_info(size_t type_size, size_t cap, size_t len, const void *data);
void mr_dyn_array_destroy(void *arr);

void *mr_dyn_array_push_alloc(void *arr, void *elem);
void mr_dyn_array_pop(void *arr);

#endif // _MR_DYN_ARRAY_H