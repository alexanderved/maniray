#ifndef _MR_MANIFOLD_H
#define _MR_MANIFOLD_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct mr_chart mr_chart;

typedef bool (*mr_chart_bounds_fn)(mr_chart *, const double *);
typedef bool (*mr_chart_period_fn)(mr_chart *, const double *);

typedef void (*mr_chart_transition_map_fn)(mr_chart *, mr_chart *, double *, const double *);

struct mr_chart {
    mr_chart_bounds_fn bounds;
    mr_chart_period_fn period;
};

typedef struct mr_manifold {
    size_t dim;

    size_t nb_charts;
    mr_chart *charts;

    mr_chart_transition_map_fn *transitions;
} mr_manifold;

#endif // _MR_MANIFOLD_H