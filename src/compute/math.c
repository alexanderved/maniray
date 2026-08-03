#include <math.h>
#include <tgmath.h>

#include "maniray/compute/math.h"

mr_float mr_norm2_2d(mr_float x, mr_float y) {
    return sqrt(x * x + y * y);
}

mr_float mr_norm2(mr_float x, mr_float y, mr_float z) {
    return sqrt(x * x + y * y + z * z);
}

mr_float mr_norm_inf(mr_float x, mr_float y, mr_float z) {
    return MR_MAX(MR_MAX(fabs(x), fabs(y)), fabs(z));
}

mr_float mr_wrap(mr_float x, mr_float min, mr_float max) {
    mr_float width = max - min;
    return fmod(fmod(x - min, width) + width, width) + min;
}

mr_float mr_atan2p(mr_float y, mr_float x) {
    return mr_wrap(atan2(y, x), 0, 2.0 * MR_PI);
}