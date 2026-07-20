#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>

#include "maniray/utils/mem_pool.h"
#include "maniray/utils/xmalloc.h"
#include "maniray/utils/misc.h"

typedef struct free_block free_block;

struct free_block {
    mr_index idx;
    mr_isize len;

    free_block *next_block;
};

typedef struct free_block_list {
    free_block *first_free_block;
    free_block *last_free_block;
} free_block_list;

struct mr_mem_pool {
    size_t capacity;

    size_t nb_types;
    size_t *block_sizes;
    size_t *offsets;
    size_t full_size;

    char *data;
    mr_uint *valid;

    free_block_list free_blocks;
};

static bool free_block_list_is_empty(free_block_list *list) {
    return !list || !list->first_free_block;
}

static free_block *find_insert_location(free_block_list *list, mr_index idx) {
    if (!list || !list->first_free_block) {
        return NULL;
    }

    free_block *b = list->first_free_block;
    while (b->next_block && b->next_block->idx < idx) {
        b = b->next_block;
    }

    if (idx < b->idx) {
        return NULL;
    }

    return b;
}

static free_block *alloc_new_block(mr_index idx, size_t len) {
    free_block *new_block = xmalloc(sizeof(free_block));

    new_block->idx = idx;
    new_block->len = len;
    new_block->next_block = NULL;

    return new_block;
}

static void free_block_list_insert(free_block_list *list, mr_index idx) {
    if (!list) {
        return;
    }

    free_block *loc = find_insert_location(list, idx);
    if (!loc) {
        free_block *n = list->first_free_block;

        if (n && n->idx == idx + 1) {
            --n->idx;
            ++n->len;
        } else if (n) {
            list->first_free_block = alloc_new_block(idx, 1);
            list->first_free_block->next_block = n;
        } else {
            list->first_free_block = list->last_free_block = alloc_new_block(idx, 1);
        }

        return;
    }

    if (loc->idx + loc->len == idx) {
        ++loc->len;
    } else if (loc->next_block && loc->next_block->idx == idx + 1) {
        --loc->next_block->idx;
        ++loc->next_block->len;
    } else {
        free_block *new_block = alloc_new_block(idx, 1);

        new_block->next_block = loc->next_block;
        loc->next_block = new_block;

        if (loc == list->last_free_block) {
            list->last_free_block = new_block;
        }

        loc = new_block;
    }

    if (loc->next_block && loc->idx + loc->len == loc->next_block->idx) {
        free_block *n = loc->next_block;

        loc->len += n->len;
        loc->next_block = n->next_block;

        free(n);
    }
}

static mr_index free_block_list_remove(free_block_list *list, mr_isize len) {
    if (!list) {
        return MR_INVALID_INDEX;
    }

    free_block *prev = NULL;
    free_block *block = list->first_free_block;
    while (block && block->len < len) {
        prev = block;
        block = block->next_block;
    }

    if (!block) {
        return MR_INVALID_INDEX;
    }

    mr_index idx = block->idx;
    block->idx += len;
    block->len -= len;

    if (block->len == 0) {
        if (prev) {
            prev->next_block = block->next_block;

            if (block == list->last_free_block) {
                list->last_free_block = prev;
            }
        } else {
            list->first_free_block = block->next_block;

            if (block == list->last_free_block) {
                list->last_free_block = list->first_free_block;
            }
        }

        free(block);
    }

    return idx;
}

mr_mem_pool *mr_mem_pool_create(size_t block_sizes[], size_t nb_types) {
    mr_mem_pool *pool = xmalloc(sizeof(mr_mem_pool));

    pool->capacity = 0;

    pool->nb_types = nb_types;
    pool->block_sizes = xmalloc(pool->nb_types * sizeof(size_t));
    memcpy(pool->block_sizes, block_sizes, pool->nb_types * sizeof(size_t));

    pool->offsets = xmalloc(pool->nb_types * sizeof(size_t));
    pool->offsets[0] = 0;
    for (size_t i = 1; i < pool->nb_types; ++i) {
        pool->offsets[i] = pool->block_sizes[i - 1] + pool->offsets[i - 1];
    }

    pool->full_size = 0;
    for (size_t i = 0; i < pool->nb_types; ++i) {
        pool->full_size += pool->block_sizes[i];
    }

    pool->data = NULL;
    pool->valid = NULL;

    pool->free_blocks = (free_block_list){ 0 };

    return pool;
}

void mr_mem_pool_destroy(mr_mem_pool *pool) {
    if (!pool) {
        return;
    }

    while (!free_block_list_is_empty(&pool->free_blocks)) {
        free_block *b = pool->free_blocks.first_free_block;
        pool->free_blocks.first_free_block = b->next_block;

        free(b);
    }

    free(pool->valid);
    free(pool->data);

    free(pool->offsets);
    free(pool->block_sizes);

    free(pool);
}

size_t mr_mem_pool_capacity(mr_mem_pool *pool) {
    return pool ? pool->capacity : 0;
}

size_t mr_mem_pool_len_bound(mr_mem_pool *pool) {
    if (!pool) {
        return 0;
    }

    free_block *last = pool->free_blocks.last_free_block;
    if (last && last->idx + last->len == (mr_isize)pool->capacity) {
        return last->idx;
    }

    return pool->capacity;
}

static bool is_element_valid(mr_mem_pool *pool, mr_index idx) {
    mr_uint chunk = pool->valid[idx / MR_INT_NB_BITS];
    mr_uint res = (chunk >> (idx % MR_INT_NB_BITS)) & 1;

    return (bool)res;
}

bool mr_mem_pool_is_block_accessible(mr_mem_pool *pool, mr_index idx) {
    return pool
        && idx != MR_INVALID_INDEX
        && idx < (mr_isize)pool->capacity
        && is_element_valid(pool, idx);
}

static void write_element(mr_mem_pool *pool, mr_index idx, va_list elem) {
    for (mr_index i = 0; i < (mr_isize)pool->nb_types; ++i) {
        void *ptr = mr_mem_pool_ptr(pool, i, idx);

        void *value = va_arg(elem, void *);
        if (value) {
            memcpy(ptr, value, pool->block_sizes[i]);
        } else {
            memset(ptr, 0, pool->block_sizes[i]);
        }
    }
}

static void change_element_validity(mr_mem_pool *pool, mr_index idx, bool value) {
    mr_uint *chunk = &pool->valid[idx / MR_INT_NB_BITS];

    if (value) {
        *chunk |= 1 << (idx % MR_INT_NB_BITS);
    } else {
        *chunk &= ~(1 << (idx % MR_INT_NB_BITS));
    }
}

static size_t closest_power_of_two(size_t n) {
    if (n == 0) {
        return 1;
    }

    --n;
    size_t nb_bits = sizeof(size_t) * CHAR_BIT;
    for (size_t s = 1; s < nb_bits; s *= 2) {
        n |= n >> s;
    }

    return ++n;
}

static size_t valid_arr_cap(size_t cap) {
    return cap / MR_INT_NB_BITS + (cap % MR_INT_NB_BITS != 0);
}

static void init_mem_pool(mr_mem_pool *pool, size_t init_cap) {
    init_cap = closest_power_of_two(init_cap);

    pool->capacity = init_cap;
    pool->data = xmalloc(pool->full_size * init_cap);

    pool->valid = xmalloc(sizeof(mr_uint) * valid_arr_cap(init_cap));
    memset(pool->valid, 0, sizeof(mr_uint) * valid_arr_cap(init_cap));
}

static void realloc_mem_pool(mr_mem_pool *pool, size_t cap_mul) {
    cap_mul = closest_power_of_two(cap_mul);

    size_t old_cap = pool->capacity;
    pool->capacity *= cap_mul;
    pool->data = xrealloc(pool->data, pool->capacity * pool->full_size);

    for (mr_index i = pool->nb_types - 1; i > 0; --i) {
        mr_index outer_idx_old = pool->offsets[i] * old_cap;
        mr_index outer_idx_new = pool->offsets[i] * pool->capacity;
        size_t size = pool->block_sizes[i] * old_cap;
        memmove(pool->data + outer_idx_new, pool->data + outer_idx_old, size);
    }

    size_t val_cap_diff = valid_arr_cap(pool->capacity) - valid_arr_cap(old_cap);
    pool->valid = xrealloc(pool->valid, sizeof(mr_uint) * valid_arr_cap(pool->capacity));
    memset(pool->valid + valid_arr_cap(old_cap), 0, sizeof(mr_uint) * val_cap_diff);

    free_block *last = pool->free_blocks.last_free_block;
    if (last && (mr_isize)old_cap == last->idx + last->len) {
        last->len += pool->capacity - old_cap;
    } else {
        free_block *block = alloc_new_block(old_cap, pool->capacity - old_cap);
        if (last) {
            pool->free_blocks.last_free_block = last->next_block = block;
        } else {
            pool->free_blocks.first_free_block = pool->free_blocks.last_free_block = block;
        }
    }
}

mr_index mr_mem_pool_alloc(mr_mem_pool *pool) {
    return mr_mem_pool_alloc_many(pool, 1);
}

mr_index mr_mem_pool_alloc_many(mr_mem_pool *pool, size_t nb_elems) {
    if (!pool || nb_elems == 0) {
        return MR_INVALID_INDEX;
    }

    mr_index idx;
    while (true) {
        if (!pool->data) {
            init_mem_pool(pool, nb_elems);
            idx = 0;

            break;
        } else if ((idx = free_block_list_remove(&pool->free_blocks, nb_elems)) != MR_INVALID_INDEX) {
            break;
        }

        size_t req = nb_elems;
        if (pool->free_blocks.last_free_block) {
            req -= pool->free_blocks.last_free_block->len;
        }

        size_t mul = req / pool->capacity + (req % pool->capacity == 0 ? 1 : 2);
        realloc_mem_pool(pool, mul);
    }

    for (size_t i = 0; i < nb_elems; ++i) {
        change_element_validity(pool, idx + i, true);
    }

    return idx;
}

mr_index mr_mem_pool_insert(mr_mem_pool *pool, ...) {
    if (!pool) {
        return MR_INVALID_INDEX;
    }

    va_list elem;
    va_start(elem, pool);

    mr_index idx = mr_mem_pool_alloc(pool);
    write_element(pool, idx, elem);

    va_end(elem);

    return idx;
}

void mr_mem_pool_remove(mr_mem_pool *pool, mr_index idx) {
    if (!mr_mem_pool_is_block_accessible(pool, idx)) {
        return;
    }

    change_element_validity(pool, idx, false);
    free_block_list_insert(&pool->free_blocks, idx);
}

void *mr_mem_pool_ptr(mr_mem_pool *pool, mr_index field, mr_index idx) {
    if (!pool || !mr_mem_pool_is_block_accessible(pool, idx)) {
        return NULL;
    }

    mr_index outer_idx = pool->offsets[field] * pool->capacity;
    mr_index inner_idx = idx * pool->block_sizes[field];

    return pool->data + outer_idx + inner_idx;
}

void *mr_mem_pool_array_ptr(mr_mem_pool *pool, mr_index field) {
    if (!pool) {
        return NULL;
    }

    return pool->data + pool->offsets[field] * pool->capacity;
}