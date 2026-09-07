#include <math.h>
#include <tgmath.h>

#include "maniray/compute/math.h"

mr_float64 mr_norm2_2d(mr_float64 x, mr_float64 y) {
    return sqrt(x * x + y * y);
}

mr_float64 mr_norm2(mr_float64 x, mr_float64 y, mr_float64 z) {
    return sqrt(x * x + y * y + z * z);
}

mr_float64 mr_norm_inf(mr_float64 x, mr_float64 y, mr_float64 z) {
    return MR_MAX(MR_MAX(fabs(x), fabs(y)), fabs(z));
}

mr_float64 mr_wrap(mr_float64 x, mr_float64 min, mr_float64 max) {
    mr_float64 width = max - min;
    return fmod(fmod(x - min, width) + width, width) + min;
}

mr_float64 mr_atan2p(mr_float64 y, mr_float64 x) {
    return mr_wrap(atan2(y, x), 0, 2.0 * MR_PI);
}