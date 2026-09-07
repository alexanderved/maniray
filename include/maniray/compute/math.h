#ifndef _MR_MATH_H
#define _MR_MATH_H

#include "maniray/utils/types.h"

#define MR_PI 3.14159265358979323846

#define MR_MIN(a, b) ((a) < (b) ? (a) : (b))
#define MR_MAX(a, b) ((a) > (b) ? (a) : (b))
#define MR_CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

#define MR_ABS(a) ((a) > 0 ? (a) : -(a))
#define MR_MOD(a, b) ((((a) % (b)) + (b)) % (b))

mr_float64 mr_norm2_2d(mr_float64 x, mr_float64 y);
mr_float64 mr_norm2(mr_float64 x, mr_float64 y, mr_float64 z);
mr_float64 mr_norm_inf(mr_float64 x, mr_float64 y, mr_float64 z);

mr_float64 mr_wrap(mr_float64 x, mr_float64 min, mr_float64 max);

// `atan2` in range [0; 2pi)
mr_float64 mr_atan2p(mr_float64 y, mr_float64 x);

#endif // _MR_MATH_H