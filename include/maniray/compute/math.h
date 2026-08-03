#ifndef _MR_MATH_H
#define _MR_MATH_H

#include "maniray/utils/types.h"

#define MR_PI 3.14159265358979323846
#define MR_MIN(a, b) ((a) < (b) ? (a) : (b))
#define MR_MAX(a, b) ((a) > (b) ? (a) : (b))

mr_float mr_norm2_2d(mr_float x, mr_float y);
mr_float mr_norm2(mr_float x, mr_float y, mr_float z);
mr_float mr_norm_inf(mr_float x, mr_float y, mr_float z);

mr_float mr_wrap(mr_float x, mr_float min, mr_float max);

// `atan2` in range [0; 2pi)
mr_float mr_atan2p(mr_float y, mr_float x);

#endif // _MR_MATH_H