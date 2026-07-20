#ifndef _MR_TYPES_H
#define _MR_TYPES_H

#include <stdint.h>
#include <stddef.h>

#define MR_INT_NB_BITS 32

typedef int32_t mr_int;
typedef uint32_t mr_uint;

typedef float mr_float;

typedef ptrdiff_t mr_isize;
typedef mr_isize mr_index;

typedef mr_uint mr_bitfield;

#endif // _MR_TYPES_H