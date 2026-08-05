#include <string.h>
#include <math.h>

#include "maniray/utils/misc.h"
#include "maniray/compute/math.h"
#include "maniray/compute/manifold.h"

#include "utils.h"

static bool transition_0_1_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return (p[1] < -0.25f && p[1] > -1.0f && mr_norm2_2d(p[0] + 2.0f, p[2]) < 1.0f)
        || (p[1] > 0.25f && p[1] < 1.0f && mr_norm2_2d(p[0] - 2.0f, p[2]) < 1.0f);
}

static int transition_0_1(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[1] > 0.0f) {
        p_out[0] -= 2.0f;
    } else {
        p_out[0] += 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_0_2_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return (p[1] < -0.25f && p[1] > -1.0f && mr_norm2_2d(p[0] - 2.0f, p[2]) < 1.0f)
        || (p[1] > 0.25f && p[1] < 1.0f && mr_norm2_2d(p[0] + 2.0f, p[2]) < 1.0f);
}

static int transition_0_2(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[1] > 0.0f) {
        p_out[0] += 2.0f;
    } else {
        p_out[0] -= 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_0_3_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    mr_float mouth_circle1[3] = { -2.0f, 0.0f, 0.0f };
    mr_float mouth_circle2[3] = { 2.0f, 0.0f, 0.0f };

    return dist_to_circle(p, mouth_circle1) < 0.75 || dist_to_circle(p, mouth_circle2) < 0.75;
}

static int transition_0_3(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float mouth_circle1[3] = { -2.0f, 0.0f, 0.0f };
    mr_float mouth_circle2[3] = { 2.0f, 0.0f, 0.0f };

    mr_float r1 = dist_to_circle(p_in, mouth_circle1);
    mr_float r2 = dist_to_circle(p_in, mouth_circle2);

    mr_float r = 0.0f;
    const mr_float *mc = NULL;
    if (r1 < 0.75f) {
        r = r1;
        mc = mouth_circle1;
    } else if (r2 < 0.75f) {
        r = r2;
        mc = mouth_circle2;
    } else {
        return MR_FAILURE;
    }

    mr_float theta = mr_atan2p(p_in[1] - mc[1], proj_radius(p_in, mc) - 1.0f);
    mr_float phi = mr_atan2p(p_in[2] - mc[2], p_in[0] - mc[0]);

    if ((r1 < 0.75f && theta > MR_PI) || r2 < 0.75f) {
        theta += 2.0f * MR_PI;
    }

    p_out[0] = r;
    p_out[1] = theta / (8.0f * MR_PI) - 0.25f;
    p_out[2] = phi / (4.0f * MR_PI) - 0.25f;

    return MR_SUCCESS;
}

static bool transition_1_0_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return p[1] < -0.25f || p[1] > 0.25f;
}

static int transition_1_0(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[1] > 0.0f) {
        p_out[0] += 2.0f;
    } else {
        p_out[0] -= 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_1_3_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return dist_to_circle(p, (mr_float[3]) { 0.0f, 0.0f, 0.0f }) < 0.75;
}

static int transition_1_3(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float mouth_circle[3] = { 0.0f, 0.0f, 0.0f };

    mr_float r = dist_to_circle(p_in, mouth_circle);
    mr_float theta = mr_atan2p(p_in[1] - mouth_circle[1], proj_radius(p_in, mouth_circle) - 1.0f) + 2.0f * MR_PI;
    mr_float phi = mr_atan2p(p_in[2] - mouth_circle[2], p_in[0] - mouth_circle[0]);

    p_out[0] = r;
    p_out[1] = theta / (8.0f * MR_PI) - 0.25f;
    p_out[2] = phi / (4.0f * MR_PI) - 0.25f;

    return MR_SUCCESS;
}

static bool transition_2_0_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return p[1] < -0.25f || p[1] > 0.25f;
}

static int transition_2_0(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));
    if (p_in[1] > 0.0f) {
        p_out[0] -= 2.0f;
    } else {
        p_out[0] += 2.0f;
    }

    return MR_SUCCESS;
}

static bool transition_2_3_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    return dist_to_circle(p, (mr_float[3]) { 0.0f, 0.0f, 0.0f }) < 0.75;
}

static int transition_2_3(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float mouth_circle[3] = { 0.0f, 0.0f, 0.0f };

    mr_float r = dist_to_circle(p_in, mouth_circle);
    mr_float theta = mr_atan2p(p_in[1] - mouth_circle[1], proj_radius(p_in, mouth_circle) - 1.0f);
    mr_float phi = mr_atan2p(p_in[2] - mouth_circle[2], p_in[0] - mouth_circle[0]);

    p_out[0] = r;
    p_out[1] = theta / (8.0f * MR_PI) - 0.25f;
    p_out[2] = phi / (4.0f * MR_PI) - 0.25f;

    return MR_SUCCESS;
}

static void chart_3_period(const mr_chart *chart, mr_float *wp, const mr_float *p) {
    MR_UNUSED(chart);

    wp[0] = p[0];
    wp[1] = mr_wrap(p[1], -0.25, 0.25);
    wp[2] = mr_wrap(p[2], -0.25, 0.25);
}

static bool transition_3_0_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    mr_float theta = (p[1] + 0.25f) * 8.0f * MR_PI;
    return theta != MR_PI && theta != 3.0f * MR_PI;
}

static int transition_3_0(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float r = p_in[0];
    mr_float theta = (p_in[1] + 0.25f) * 8.0f * MR_PI;
    mr_float phi = (p_in[2] + 0.25f) * 4.0f * MR_PI;

    mr_float x_offset = 0.0f;
    if ((theta >= 0.0f && theta < 1.0f * MR_PI) || (theta > 3.0f * MR_PI && theta <= 4.0f * MR_PI)) {
        x_offset = -2.0f;
    } else if (theta > 1.0f * MR_PI && theta < 3.0f * MR_PI) {
        x_offset = 2.0f;
    } else {
        return MR_FAILURE;
    }

    p_out[0] = (1.0f + r * cosf(theta)) * cosf(phi) + x_offset;
    p_out[1] = r * sinf(theta);
    p_out[2] = (1.0f + r * cosf(theta)) * sinf(phi);

    return MR_SUCCESS;
}

static bool transition_3_1_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    mr_float theta = (p[1] + 0.25f) * 8.0f * MR_PI;
    return theta > 2.5f * MR_PI && theta < 3.5f * MR_PI;
}

static int transition_3_1(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float r = p_in[0];
    mr_float theta = (p_in[1] + 0.25f) * 8.0f * MR_PI;
    mr_float phi = (p_in[2] + 0.25f) * 4.0f * MR_PI;

    p_out[0] = (1.0f + r * cosf(theta)) * cosf(phi);
    p_out[1] = r * sinf(theta);
    p_out[2] = (1.0f + r * cosf(theta)) * sinf(phi);

    return MR_SUCCESS;
}

static bool transition_3_2_domain(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);

    mr_float theta = (p[1] + 0.25f) * 8.0f * MR_PI;
    return theta > 0.5f * MR_PI && theta < 1.5f * MR_PI;
}

static int transition_3_2(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);

    mr_float r = p_in[0];
    mr_float theta = (p_in[1] + 0.25f) * 8.0f * MR_PI;
    mr_float phi = (p_in[2] + 0.25f) * 4.0f * MR_PI;

    p_out[0] = (1.0f + r * cosf(theta)) * cosf(phi);
    p_out[1] = r * sinf(theta);
    p_out[2] = (1.0f + r * cosf(theta)) * sinf(phi);

    return MR_SUCCESS;
}