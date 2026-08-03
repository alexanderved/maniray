#include <string.h>

#include "maniray/compute/manifold.h"
#include "maniray/utils/misc.h"
#include "maniray/utils/xmalloc.h"

mr_float mr_chart_euclidean_metric(const mr_chart *chart, const mr_float *p, size_t i, size_t j) {
    MR_UNUSED(chart);
    MR_UNUSED(p);

    return i == j ? 1.0f : 0.0f;
}

void mr_chart_default_dtor(mr_chart_desc *chart) {
    free(chart->userdata);
}

mr_float mr_inner_product(const mr_manifold *manifold, size_t chart_idx, const mr_float *p, const mr_float *v1, const mr_float *v2) {
    mr_float res = 0.0f;
    for (size_t i = 0; i < manifold->dim; ++i) {
        for (size_t j = 0; j < manifold->dim; ++j) {
            res += mr_manifold_metric(manifold, chart_idx, p, i, j) * v1[i] * v2[j];
        }
    }

    return res;
}

static bool mr_transition_domain_self(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);
    MR_UNUSED(p);

    return true;
}

static int mr_transition_map_self(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    memcpy(p_out, p_in, t->manifold->dim * sizeof(mr_float));

    return MR_SUCCESS;
}

mr_transition_desc mr_transition_desc_create_self() {
    return (mr_transition_desc) { .domain = mr_transition_domain_self, .fn = mr_transition_map_self };
}

static bool mr_transition_domain_empty(const mr_transition *t, const mr_float *p) {
    MR_UNUSED(t);
    MR_UNUSED(p);

    return false;
}

static int mr_transition_map_empty(const mr_transition *t, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(t);
    MR_UNUSED(p_out);
    MR_UNUSED(p_in);

    return MR_FAILURE;
}

mr_transition_desc mr_transition_desc_create_empty() {
    return (mr_transition_desc) { .domain = mr_transition_domain_empty, .fn = mr_transition_map_empty };
}

mr_manifold *mr_manifold_create(size_t dim, size_t nb_charts, const mr_chart_desc *charts, const mr_transition_desc *transitions) {
    mr_manifold *manifold = xmalloc(sizeof(mr_manifold));

    manifold->dim = dim;
    manifold->nb_charts = nb_charts;

    manifold->charts = xcalloc(nb_charts, sizeof(mr_chart_desc));
    for (size_t i = 0; i < nb_charts; ++i) {
        manifold->charts[i] = (mr_chart_desc) {
            .bounds = charts[i].bounds,
            .period = charts[i].period,
            .metric = charts[i].metric ? charts[i].metric : mr_chart_euclidean_metric,
            .userdata = charts[i].userdata,
            .dtor = charts[i].dtor ? charts[i].dtor : mr_chart_default_dtor,
        };
    }

    manifold->transitions = xcalloc(nb_charts * nb_charts, sizeof(mr_transition_desc));
    memcpy(manifold->transitions, transitions, nb_charts * nb_charts * sizeof(mr_transition_desc));

    return manifold;
}

void mr_manifold_destroy(mr_manifold *manifold) {
    free(manifold->transitions);
    
    mr_chart_desc *chart = NULL;
    for (size_t i = 0; i < manifold->nb_charts; ++i) {
        chart = &manifold->charts[i];
        chart->dtor(chart);
    }
    free(manifold->charts);

    free(manifold);
}

static mr_chart mr_manifold_get_chart(const mr_manifold *manifold, size_t chart_idx) {
    mr_chart_desc *desc = &manifold->charts[chart_idx];
    return (mr_chart) {
        .manifold = manifold,
        .bounds = desc->bounds,
        .period = desc->period,
        .metric = desc->metric,
        .userdata = desc->userdata,
    };
}

#define MAX_STACK_ARR_SIZE 4

#define TEMP_ALLOC(ptr, type, size) \
    type _temp_arr_ ## ptr[MAX_STACK_ARR_SIZE] = { 0 }; \
    type *ptr = size <= MAX_STACK_ARR_SIZE ? _temp_arr_ ## ptr : xmalloc(size * sizeof(type))

#define TEMP_FREE(ptr, size) if (size > MAX_STACK_ARR_SIZE) free(ptr)

bool mr_manifold_is_in_bounds(const mr_manifold *manifold, size_t chart_idx, const mr_float *p) {
    mr_chart chart = mr_manifold_get_chart(manifold, chart_idx);

    TEMP_ALLOC(wp, mr_float, manifold->dim);
    mr_manifold_periodic_wrap(manifold, chart_idx, wp, p);

    bool res = chart.bounds(&chart, wp);
    TEMP_FREE(wp, manifold->dim);

    return res;
}

void mr_manifold_periodic_wrap(const mr_manifold *manifold, size_t chart_idx, mr_float *wp, const mr_float *p) {
    mr_chart chart = mr_manifold_get_chart(manifold, chart_idx);
    if (chart.period) {
        chart.period(&chart, wp, p);
    } else {
        memcpy(wp, p, manifold->dim * sizeof(mr_float));
    }
}

mr_float mr_manifold_metric(const mr_manifold *manifold, size_t chart_idx, const mr_float *p, size_t i, size_t j) {
    mr_chart chart = mr_manifold_get_chart(manifold, chart_idx);

    TEMP_ALLOC(wp, mr_float, manifold->dim);
    mr_manifold_periodic_wrap(manifold, chart_idx, wp, p);

    mr_float res = chart.metric(&chart, wp, i, j);
    TEMP_FREE(wp, manifold->dim);

    return res;
}

int mr_manifold_transition(const mr_manifold *manifold, size_t i, size_t j, mr_float *p_out, const mr_float *p_in) {
    mr_chart src = mr_manifold_get_chart(manifold, i);
    mr_chart dst = mr_manifold_get_chart(manifold, j);

    mr_transition_desc *desc = &manifold->transitions[i * manifold->nb_charts + j];
    mr_transition t = {
        .manifold = manifold,
        .src = &src,
        .dst = &dst,
        .domain = desc->domain,
        .fn = desc->fn,
    };

    TEMP_ALLOC(wp, mr_float, manifold->dim);
    mr_manifold_periodic_wrap(manifold, i, wp, p_in);

    int res = MR_FAILURE;
    if (t.domain(&t, wp) && t.src->bounds(t.src, wp)) {
        res = t.fn(&t, p_out, wp);
        mr_manifold_periodic_wrap(manifold, j, p_out, p_out);
    }

    TEMP_FREE(wp, manifold->dim);

    return res;
}