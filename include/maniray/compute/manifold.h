#ifndef _MR_MANIFOLD_H
#define _MR_MANIFOLD_H

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

#include "maniray/utils/types.h"

typedef struct mr_chart mr_chart;
typedef struct mr_chart_desc mr_chart_desc;

typedef struct mr_manifold mr_manifold;
typedef struct mr_transition mr_transition;

typedef bool (*mr_chart_bounds_fn)(const mr_chart *chart, const mr_float *p);
typedef void (*mr_chart_period_fn)(const mr_chart *chart, mr_float *wp, const mr_float *p);

typedef mr_float (*mr_chart_tensor2d_fn)(const mr_chart *chart, const mr_float *p, size_t i, size_t j);

struct mr_chart_desc {
    mr_chart_bounds_fn bounds;
    mr_chart_period_fn period;

    mr_chart_tensor2d_fn metric;
    mr_chart_tensor2d_fn inv_metric;
};

struct mr_chart {
    const mr_manifold *manifold;
    size_t idx;

    mr_chart_bounds_fn bounds;
    mr_chart_period_fn period;

    mr_chart_tensor2d_fn metric;
    mr_chart_tensor2d_fn inv_metric;
};

mr_float mr_chart_euclidean_metric(const mr_chart *chart, const mr_float *p, size_t i, size_t j);
mr_float mr_inner_product(const mr_manifold *manifold, size_t chart_idx, const mr_float *p, const mr_float *v1, const mr_float *v2);

typedef bool (*mr_transition_domain_fn)(const mr_transition *transition, const mr_float *p);
typedef int (*mr_transition_map_fn)(const mr_transition *transition, mr_float *p_out, const mr_float *p_in);

typedef struct mr_transition_desc {
    mr_transition_domain_fn domain;
    mr_transition_map_fn fn;
} mr_transition_desc;

mr_transition_desc mr_transition_desc_create_self();
mr_transition_desc mr_transition_desc_create_empty();

struct mr_transition {
    const mr_manifold *manifold;

    const mr_chart *src;
    const mr_chart *dst;

    mr_transition_domain_fn domain;
    mr_transition_map_fn fn;
};

typedef void (*mr_manifold_userdata_dtor_fn)(void *userdata);

struct mr_manifold {
    size_t dim;

    size_t nb_charts;
    mr_chart_desc *charts;

    mr_transition_desc *transitions;

    void *userdata;
    mr_manifold_userdata_dtor_fn dtor;
};

mr_manifold *mr_manifold_create(
    size_t dim,
    size_t nb_charts,
    const mr_chart_desc *charts,
    const mr_transition_desc *transitions
);
mr_manifold *mr_manifold_create_with_userdata(
    size_t dim,
    size_t nb_charts,
    const mr_chart_desc *charts,
    const mr_transition_desc *transitions,
    void *userdata,
    mr_manifold_userdata_dtor_fn dtor
);

void mr_manifold_destroy(mr_manifold *manifold);

void mr_manifold_userdata_default_dtor(void *userdata);
void mr_manifold_userdata_noop_dtor(void *userdata);

bool mr_manifold_is_in_bounds(const mr_manifold *manifold, size_t chart_idx, const mr_float *p);
void mr_manifold_periodic_wrap(const mr_manifold *manifold, size_t chart_idx, mr_float *wp, const mr_float *p);

mr_float mr_manifold_metric(const mr_manifold *manifold, size_t chart_idx, const mr_float *p, size_t i, size_t j);
mr_float mr_manifold_inv_metric(const mr_manifold *manifold, size_t chart_idx, const mr_float *p, size_t i, size_t j);

int mr_manifold_transition(const mr_manifold *manifold, size_t i, size_t j, mr_float *p_out, const mr_float *p_in);

#endif // _MR_MANIFOLD_H