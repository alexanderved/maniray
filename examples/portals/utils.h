#ifndef _PORTALS_UTILS_H
#define _PORTALS_UTILS_H

#include "maniray/utils/types.h"
#include "maniray/compute/math.h"

static mr_float proj_radius(const mr_float *p, const mr_float *c) {
    return mr_norm2_2d(p[0] - c[0], p[2] - c[2]);
}

static mr_float dist_to_circle(const mr_float *p, const mr_float *c) {
    mr_float proj = proj_radius(p, c) - 1.0;
    mr_float n = p[1] - c[1];

    return mr_norm2_2d(proj, n);
}

#endif // _PORTALS_UTILS_H