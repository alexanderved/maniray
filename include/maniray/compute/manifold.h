#ifndef _MR_MANIFOLD_H
#define _MR_MANIFOLD_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "maniray/utils/types.h"

typedef struct mr_chart mr_chart;
typedef struct mr_manifold mr_manifold;

typedef bool (*mr_chart_bounds_fn)(const mr_manifold *, const mr_chart *, const mr_float *);
typedef mr_float (*mr_chart_metric_fn)(const mr_manifold *, const mr_chart *, const mr_float *, size_t, size_t);

typedef void (*mr_chart_transition_map_fn)(const mr_manifold *, const mr_chart *, const mr_chart *, mr_float *, const mr_float *);

mr_float mr_chart_euclidean_metric(const mr_manifold *manifold, const mr_chart *chart, const mr_float *p, size_t i, size_t j);
void mr_chart_transition_dummy(const mr_manifold *manifold, const mr_chart *c1, const mr_chart *c2, mr_float *p_out, const mr_float *p_in);
mr_float mr_chart_inner_product(const mr_manifold *manifold, const mr_chart *chart, const mr_float *p, const mr_float *v1, const mr_float *v2);

struct mr_chart {
    mr_chart_bounds_fn bounds;
    mr_chart_metric_fn metric;
};

void mr_chart_init(mr_chart *chart, mr_chart_bounds_fn bounds, mr_chart_metric_fn metric);

struct mr_manifold {
    size_t dim;

    size_t nb_charts;
    mr_chart *charts;

    mr_chart_transition_map_fn *transitions;
};

mr_manifold *mr_manifold_create(size_t dim, size_t nb_charts);
void mr_manifold_destroy(mr_manifold *manifold);

mr_chart *mr_manifold_get_chart(mr_manifold *manifold, size_t i);
void mr_manifold_set_chart(mr_manifold *manifold, size_t i, const mr_chart *chart);

mr_chart_transition_map_fn mr_manifold_get_transition(mr_manifold *manifold, size_t i, size_t j);
void mr_manifold_set_transition(mr_manifold *manifold, size_t i, size_t j, mr_chart_transition_map_fn map);

#endif // _MR_MANIFOLD_H