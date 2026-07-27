#include <string.h>

#include "maniray/compute/manifold.h"
#include "maniray/utils/misc.h"
#include "maniray/utils/xmalloc.h"

mr_float mr_chart_euclidean_metric(const mr_manifold *manifold, const mr_chart *chart, const mr_float *p, size_t i, size_t j) {
    MR_UNUSED(manifold);
    MR_UNUSED(chart);
    MR_UNUSED(p);

    return i == j ? 1.0 : 0.0;
}

void mr_chart_transition_dummy(const mr_manifold *manifold, const mr_chart *c1, const mr_chart *c2, mr_float *p_out, const mr_float *p_in) {
    MR_UNUSED(c1);
    MR_UNUSED(c2);

    memcpy(p_out, p_in, sizeof(mr_float) * manifold->dim);
}

mr_float mr_chart_inner_product(const mr_manifold *manifold, const mr_chart *chart, const mr_float *p, const mr_float *v1, const mr_float *v2) {
    mr_float res = 0.0f;
    for (size_t i = 0; i < manifold->dim; ++i) {
        for (size_t j = 0; j < manifold->dim; ++j) {
            res += chart->metric(manifold, chart, p, i, j) * v1[i] * v2[j];
        }
    }

    return res;
}

void mr_chart_init(mr_chart *chart, mr_chart_bounds_fn bounds, mr_chart_metric_fn metric) {
    chart->bounds = bounds;
    chart->metric = metric;
}

mr_manifold *mr_manifold_create(size_t dim, size_t nb_charts) {
    mr_manifold *manifold = xmalloc(sizeof(mr_manifold));

    manifold->dim = dim;
    manifold->nb_charts = nb_charts;

    manifold->charts = xcalloc(nb_charts, sizeof(mr_chart));
    manifold->transitions = xcalloc(nb_charts * nb_charts, sizeof(mr_chart_transition_map_fn));

    return manifold;
}

void mr_manifold_destroy(mr_manifold *manifold) {
    free(manifold->transitions);
    free(manifold->charts);
    free(manifold);
}

mr_chart *mr_manifold_get_chart(mr_manifold *manifold, size_t i) {
    return &manifold->charts[i];
}

void mr_manifold_set_chart(mr_manifold *manifold, size_t i, const mr_chart *chart) {
    manifold->charts[i] = *chart;
}

mr_chart_transition_map_fn mr_manifold_get_transition(mr_manifold *manifold, size_t i, size_t j) {
    mr_chart_transition_map_fn t = *(manifold->transitions + i * manifold->nb_charts + j);

    return t ? t : mr_chart_transition_dummy;
}

void mr_manifold_set_transition(mr_manifold *manifold, size_t i, size_t j, mr_chart_transition_map_fn map) {
    *(manifold->transitions + i * manifold->nb_charts + j) = map;
}