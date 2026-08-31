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

static mr_float chart_3_metric(const mr_chart *chart, const mr_float *p, size_t i, size_t j) {
    MR_UNUSED(chart);

    if (i != j) {
        return 0.0f;
    }

    mr_float res_sqrt = 0.0f;
    switch (i) {
        case 0:
            res_sqrt = 1.0f;
            break;

        case 1: ;
            res_sqrt = 8.0f * MR_PI * p[0];
            break;

        case 2: ;
            res_sqrt = 4.0f * MR_PI * (1 + p[0] * cosf(8.0f * MR_PI * p[1]));
            break;
    }

    return res_sqrt * res_sqrt;
}

static mr_float chart_3_inverse_metric(const mr_chart *chart, const mr_float *p, size_t i, size_t j) {
    MR_UNUSED(chart);

    return i == j ? 1.0f / chart_3_metric(chart, p, i, j) : 0.0f;
}