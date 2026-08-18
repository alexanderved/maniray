#include "maniray/compute/manifold3d.h"
#include "maniray/utils/types.h"

mr_float mr_manifold3d_metric_minor(const mr_manifold *manifold, size_t chart_idx, const mr_float p[MR_NB_AXES], size_t exc_i, size_t exc_j) {
    if (!manifold || !p) {
        return 0.0f;
    }

    size_t curr_idx = 0;
    mr_float submatrix[4] = { 0.0f };
    for (size_t i = 0; i < MR_NB_AXES; ++i) {
        if (i == exc_i) {
            continue;
        }

        for (size_t j = 0; j < MR_NB_AXES; ++j) {
            if (j == exc_j) {
                continue;
            }

            submatrix[curr_idx++] = mr_manifold_metric(manifold, chart_idx, p, i, j);
        }
    }

    return submatrix[0] * submatrix[3] - submatrix[1] * submatrix[2];
}

mr_float mr_manifold3d_metric_determinant(const mr_manifold *manifold, size_t chart_idx, const mr_float p[MR_NB_AXES]) {
    if (!manifold || !p) {
        return 0.0f;
    }

    mr_float det = 0.0f;
    for (size_t j = 0; j < MR_NB_AXES; ++j) {
        det += (1 - j % 2 * 2)
            * mr_manifold_metric(manifold, chart_idx, p, 0, j)
            * mr_manifold3d_metric_minor(manifold, chart_idx, p, 0, j);
    }

    return det;
}

mr_float mr_manifold3d_chart_auto_inv_metric(const mr_chart *chart, const mr_float p[MR_NB_AXES], size_t i, size_t j) {
    if (!chart || !p) {
        return 0.0f;
    }

    mr_float sign = 1 - (i + j) % 2 * 2;
    mr_float minor = mr_manifold3d_metric_minor(chart->manifold, chart->idx, p, j, i);
    mr_float det = mr_manifold3d_metric_determinant(chart->manifold, chart->idx, p);

    return sign * minor / det;
}