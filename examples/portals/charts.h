#include "maniray/utils/misc.h"
#include "maniray/compute/math.h"
#include "maniray/compute/manifold.h"

static bool chart_0_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);

    return p[1] > 0.25f || p[1] < -0.25f || (mr_norm2_2d(p[0] + 2.0f, p[2]) > 1.0f && mr_norm2_2d(p[0] - 2.0f, p[2]) > 1.0f);
}

static bool chart_1_2_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);

    return p[1] < 1.0f && p[1] > -1.0f && mr_norm2_2d(p[0], p[2]) < 1.0f;
}

static bool chart_3_bounds(const mr_chart *chart, const mr_float *p) {
    MR_UNUSED(chart);

    return p[0] > 0.0 && p[0] < 0.75;
}