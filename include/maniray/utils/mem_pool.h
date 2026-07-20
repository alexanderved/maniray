#ifndef _MR_MEM_POOL_H
#define _MR_MEM_POOL_H

#include <stdlib.h>
#include <stdbool.h>

#include "maniray/utils/types.h"

#define MR_INVALID_INDEX -1

typedef struct mr_mem_pool mr_mem_pool;

mr_mem_pool *mr_mem_pool_create(size_t block_sizes[], size_t nb_types);
void mr_mem_pool_destroy(mr_mem_pool *pool);

size_t mr_mem_pool_capacity(mr_mem_pool *pool);
size_t mr_mem_pool_len_bound(mr_mem_pool *pool);

bool mr_mem_pool_is_block_accessible(mr_mem_pool *pool, mr_index idx);

mr_index mr_mem_pool_alloc(mr_mem_pool *pool);
mr_index mr_mem_pool_alloc_many(mr_mem_pool *pool, size_t nb_elems);

mr_index mr_mem_pool_insert(mr_mem_pool *pool, ...);
void mr_mem_pool_remove(mr_mem_pool *pool, mr_index idx);

void *mr_mem_pool_array_ptr(mr_mem_pool *pool, mr_index field);
void *mr_mem_pool_ptr(mr_mem_pool *pool, mr_index field, mr_index idx);

#endif // _MR_MEM_POOL_H