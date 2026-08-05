#include "charts.h"
#include "transitions.h"

static mr_manifold *setup_manifold() {
#define NB_CHARTS 4
    mr_chart_desc charts[NB_CHARTS] = {
        { .bounds = chart_0_bounds },
        { .bounds = chart_1_2_bounds },
        { .bounds = chart_1_2_bounds },
        { .bounds = chart_3_bounds, .period = chart_3_period }, // TODO: Add metric
    };

    mr_transition_desc transitions[NB_CHARTS * NB_CHARTS] = {
        mr_transition_desc_create_self(),
        { transition_0_1_domain, transition_0_1 },
        { transition_0_2_domain, transition_0_2 },
        { transition_0_3_domain, transition_0_3 },

        { transition_1_0_domain, transition_1_0 },
        mr_transition_desc_create_self(),
        mr_transition_desc_create_empty(),
        { transition_1_3_domain, transition_1_3 },

        { transition_2_0_domain, transition_2_0 },
        mr_transition_desc_create_empty(),
        mr_transition_desc_create_self(),
        { transition_2_3_domain, transition_2_3 },

        { transition_3_0_domain, transition_3_0 },
        { transition_3_1_domain, transition_3_1 },
        { transition_3_2_domain, transition_3_2 },
        mr_transition_desc_create_self(),
    };

    return mr_manifold_create(3, NB_CHARTS, charts, transitions);
}