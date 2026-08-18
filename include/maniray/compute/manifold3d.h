#ifndef _MR_MANIFOLD_3D_H
#define _MR_MANIFOLD_3D_H

#include "maniray/compute/manifold.h"
#include "maniray/compute/navigation.h"

mr_float mr_manifold3d_metric_minor(const mr_manifold *manifold, size_t chart_idx, const mr_float p[MR_NB_AXES], size_t exc_i, size_t exc_j);
mr_float mr_manifold3d_metric_determinant(const mr_manifold *manifold, size_t chart_idx, const mr_float p[MR_NB_AXES]);

mr_float mr_manifold3d_chart_auto_inv_metric(const mr_chart *chart, const mr_float p[MR_NB_AXES], size_t i, size_t j);

#endif // _MR_MANIFOLD_3D_H