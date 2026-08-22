#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "maniray/compute/math.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/utils/misc.h"

// Q for triQuadratic interpolation
#define Q_STENCIL_DIM 3
#define Q_STENCIL_DIM_SQR (Q_STENCIL_DIM * Q_STENCIL_DIM)

static mr_int extract_coord(mr_int local_idx, mr_axis axis) {
    return local_idx >> (axis * 2) & 0x3;
}

static mr_int pack_local_idx(mr_int coords[MR_NB_AXES]) {
    return (coords[2] << 4) | (coords[1] << 2) | coords[0];
}

static mr_float lagrange_coef(mr_float x, mr_float stencil[Q_STENCIL_DIM], mr_int idx) {
    mr_float res = 1.0f;
    for (mr_int i = 0; i < Q_STENCIL_DIM; ++i) {
        if (i != idx) {
            res *= (x - stencil[i]) / (stencil[idx] - stencil[i]);
        }
    }

    return res;
}

static mr_float lagrange_coef_with_idx(mr_float x, mr_int first_idx, mr_int idx) {
    mr_float stencil[] = {
        (mr_float)first_idx,
        (mr_float)(first_idx + 1),
        (mr_float)(first_idx + 2),
    };

    return lagrange_coef(x, stencil, idx - first_idx);
}

static mr_float calculate_coef(
    const mr_float p[MR_NB_AXES],
    mr_octree_cell *cell,
    mr_int first_local_idx[MR_NB_AXES],
    mr_int local_idx[MR_NB_AXES],
    mr_int host_local_idx[MR_NB_AXES]
) {
    mr_float local_coords[] = {
        (p[MR_AXIS_X] - cell->x) / cell->dim + host_local_idx[MR_AXIS_X],
        (p[MR_AXIS_Y] - cell->y) / cell->dim + host_local_idx[MR_AXIS_Y],
        (p[MR_AXIS_Z] - cell->z) / cell->dim + host_local_idx[MR_AXIS_Z],
    };

    mr_float coef = 1.0f;
    for (mr_int i = 0; i < Q_STENCIL_DIM; ++i) {
        coef *= lagrange_coef_with_idx(local_coords[i], first_local_idx[i], local_idx[i]);
    }

    return coef;
}

static mr_int get_first_cell_local_coord(mr_int local_coord) {
    return local_coord / 2;
}

int mr_fvm_perform_interpolation(mr_ocforest *forest, mr_index octree_idx, const mr_float p[MR_NB_AXES], mr_fvm_interpolation_cb interp) {
    assert(forest);
    assert(p);
    assert(interp.fn);
    assert((size_t)octree_idx < forest->nb_roots);

    mr_int host_cell_idx = mr_octree_locate_point_in_cell(forest, octree_idx, p);
    if (host_cell_idx == MR_INVALID_INDEX) {
        return MR_FAILURE;
    }

    mr_octree_cell *host_cell = mr_ocforest_get_cell(forest, host_cell_idx);
    mr_octree_node *host_node = mr_ocforest_get_node(forest, host_cell->parent);
    if (!(host_node->flags & MR_OCTREE_NODE_FLAG_ACTIVE)) {
        return MR_FAILURE;
    }

    mr_int host_local_idx = host_cell_idx - host_node->first_child;
    mr_int host_local_coords[] = {
        extract_coord(host_local_idx, MR_AXIS_X),
        extract_coord(host_local_idx, MR_AXIS_Y),
        extract_coord(host_local_idx, MR_AXIS_Z),
    };

    mr_int first_local_coords[] = {
        get_first_cell_local_coord(host_local_coords[0]),
        get_first_cell_local_coord(host_local_coords[1]),
        get_first_cell_local_coord(host_local_coords[2]),
    };

    for (mr_int i = 0; i < Q_STENCIL_DIM; ++i) {
        for (mr_int j = 0; j < Q_STENCIL_DIM; ++j) {
            for (mr_int k = 0; k < Q_STENCIL_DIM; ++k) {
                mr_int local_coords[] = {
                    first_local_coords[0] + i,
                    first_local_coords[1] + j,
                    first_local_coords[2] + k,
                };

                mr_int cell_idx = host_node->first_child + pack_local_idx(local_coords);
                mr_float coef = calculate_coef(p, host_cell, first_local_coords, local_coords, host_local_coords);

                if (interp.fn(forest, cell_idx, coef, interp.userdata) != MR_SUCCESS) {
                    return MR_FAILURE;
                }
            }
        }
    }

    return MR_SUCCESS;
}

static mr_float get_cell_center_axis(const mr_octree_cell *cell, mr_axis axis) {
    switch (axis) {
        case MR_AXIS_X:
            return cell->x;
        case MR_AXIS_Y:
            return cell->y;
        case MR_AXIS_Z:
            return cell->z;
        default:
            assert(false);
            return 0.0f;
    }
}

static mr_direction calculate_ghost_cell_direction(mr_octree_cell *coarse_cell, mr_octree_cell *fine_cell) {
    for (mr_axis axis = MR_AXIS_X; axis <= MR_AXIS_Z; ++axis) {
        mr_float diff = get_cell_center_axis(fine_cell, axis) - get_cell_center_axis(coarse_cell, axis);
        if (MR_ABS(diff) > coarse_cell->dim / 2.0f) {
            mr_sign sign = diff < 0.0 ? MR_SIGN_MINUS : MR_SIGN_PLUS;

            return mr_direction_create(axis, sign);
        }
    }

    assert(false);
    return MR_DIRECTION_MI_X;
}

static mr_int get_tangent_stencil_first_local_idx(
    mr_ocforest *forest,
    mr_int coarse_cell_idx,
    mr_octree_cell *cell,
    mr_direction dir
) {
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_int local_idx = coarse_cell_idx - node->first_child;

    mr_int local_coords[] = {
        extract_coord(local_idx, MR_AXIS_X),
        extract_coord(local_idx, MR_AXIS_Y),
        extract_coord(local_idx, MR_AXIS_Z),
    };

    mr_axis normal_axis = mr_direction_get_axis(dir);
    mr_int first_local_coords[MR_NB_AXES] = { 0 };

    for (mr_axis axis = MR_AXIS_X; axis <= MR_AXIS_Z; ++axis) {
        first_local_coords[axis] = (axis == normal_axis) ? local_coords[axis]
                                                         : get_first_cell_local_coord(local_coords[axis]);
    }

    return pack_local_idx(first_local_coords);
}

static void calculate_tangent_data_point(
    mr_octree_cell *coarse_cell,
    mr_octree_cell *fine_cell,
    mr_direction dir,
    mr_float tangent_data_point[MR_NB_AXES]
) {
    mr_axis normal_axis = mr_direction_get_axis(dir);

    for (mr_axis axis = MR_AXIS_X; axis <= MR_AXIS_Z; ++axis) {
        tangent_data_point[axis] = (axis == normal_axis) ? get_cell_center_axis(coarse_cell, axis)
                                                         : get_cell_center_axis(fine_cell, axis);
    }
}

static mr_float calculate_tangent_coef(
    const mr_float p[2],
    mr_int first_local_coords[2],
    mr_int local_coords[2]
) {
    return lagrange_coef_with_idx(p[0], first_local_coords[0], local_coords[0])
        * lagrange_coef_with_idx(p[1], first_local_coords[1], local_coords[1]);
}

static void get_normal_stencil_indices(mr_ocforest *forest, mr_int fine_cell_idx, mr_direction dir, mr_int stencil[2]) {
    mr_octree_cell_neighbor n = mr_octree_find_face_neighbor_cells(forest, fine_cell_idx, dir);
    assert(n.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);

    stencil[0] = n.neighbor_idx;
    stencil[1] = fine_cell_idx;
}

static void calculate_normal_stencil(
    mr_ocforest *forest,
    mr_int stencil_cell_indices[2],
    mr_float tangent_data_point[MR_NB_AXES],
    mr_direction dir,
    mr_float stencil[Q_STENCIL_DIM]
) {
    mr_axis normal_axis = mr_direction_get_axis(dir);

    mr_octree_cell *cells[] = {
        mr_ocforest_get_cell(forest, stencil_cell_indices[0]),
        mr_ocforest_get_cell(forest, stencil_cell_indices[1]),
    };

    stencil[0] = get_cell_center_axis(cells[0], normal_axis);
    stencil[1] = get_cell_center_axis(cells[1], normal_axis);
    stencil[2] = tangent_data_point[normal_axis];
}

static void calculate_normal_coefs(mr_float stencil[Q_STENCIL_DIM], mr_float coefs[Q_STENCIL_DIM]) {
    // The distance from the ghost cell to the finer cell center is cell->dim / 2,
    // to the coarser cell center is cell->dim / 4. Therefore, the ghost cell divides
    // the interval between points in a 2:1 ratio.
    mr_float ghost_cell_coord = (stencil[1] + 2.0f * stencil[2]) / 3.0f; 

    for (size_t i = 0; i < Q_STENCIL_DIM; ++i) {
        coefs[i] = lagrange_coef(ghost_cell_coord, stencil, i);
    }
}

static int interpolate_ghost_cell_tangent(
    mr_ocforest *forest,
    mr_int coarse_cell_idx,
    mr_int tangent_first_local_idx,
    mr_float tangent_data_point[MR_NB_AXES],
    mr_direction dir,
    mr_float mul,
    mr_fvm_interpolation_cb interp
) {
    mr_octree_cell *cell = mr_ocforest_get_cell(forest, coarse_cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);

    mr_axis normal_axis = mr_direction_get_axis(dir);
    mr_axis tangent_axes[2] = {
        MR_MOD((mr_int)normal_axis - 1, MR_NB_AXES),
        MR_MOD((mr_int)normal_axis + 1, MR_NB_AXES),
    };

    mr_int cell_local_idx = coarse_cell_idx - node->first_child;
    mr_int cell_tangent_coords[] = {
        extract_coord(cell_local_idx, tangent_axes[0]),
        extract_coord(cell_local_idx, tangent_axes[1]),
    };

    mr_int first_coord_normal = extract_coord(tangent_first_local_idx, normal_axis);
    mr_int first_coords_tangent[] = {
        extract_coord(tangent_first_local_idx, tangent_axes[0]),
        extract_coord(tangent_first_local_idx, tangent_axes[1]),
    };

    mr_float tangent_data_point_coords[] = {
        (tangent_data_point[tangent_axes[0]] - get_cell_center_axis(cell, tangent_axes[0])) / cell->dim + cell_tangent_coords[0],
        (tangent_data_point[tangent_axes[1]] - get_cell_center_axis(cell, tangent_axes[1])) / cell->dim + cell_tangent_coords[1],
    };

    for (mr_int i = 0; i < MR_NB_AXES; ++i) {
        for (mr_int j = 0; j < MR_NB_AXES; ++j) {
            mr_int tangent_coords[] = {
                first_coords_tangent[0] + i,
                first_coords_tangent[1] + j,
            };

            mr_int local_coords[3] = { 0 };
            local_coords[normal_axis] = first_coord_normal;
            local_coords[tangent_axes[0]] = tangent_coords[0];
            local_coords[tangent_axes[1]] = tangent_coords[1];

            mr_int cell_idx = node->first_child + pack_local_idx(local_coords);
            mr_float coef = mul * calculate_tangent_coef(tangent_data_point_coords, first_coords_tangent, tangent_coords);

            if (interp.fn(forest, cell_idx, coef, interp.userdata) != MR_SUCCESS) {
                return MR_FAILURE;
            }
        }
    }

    return MR_SUCCESS;
}

static int interpolate_ghost_cell(
    mr_ocforest *forest,
    mr_int coarse_cell_idx,
    mr_int tangent_first_local_idx,
    mr_int normal_stencil_cell_indices[2],
    mr_float tangent_data_point[MR_NB_AXES],
    mr_direction dir,
    mr_fvm_interpolation_cb interp
) {
    mr_float normal_stencil[Q_STENCIL_DIM] = { 0.0f };
    calculate_normal_stencil(forest, normal_stencil_cell_indices, tangent_data_point, dir, normal_stencil);

    mr_float normal_coefs[Q_STENCIL_DIM] = { 0.0f };
    calculate_normal_coefs(normal_stencil, normal_coefs);

    for (size_t i = 0; i < 2; ++i) {
        mr_int cell_idx = normal_stencil_cell_indices[i];
        mr_float coef = normal_coefs[i];

        if (interp.fn(forest, cell_idx, coef, interp.userdata) != MR_SUCCESS) {
            return MR_FAILURE;
        }
    }

    return interpolate_ghost_cell_tangent(
        forest,
        coarse_cell_idx,
        tangent_first_local_idx,
        tangent_data_point,
        dir,
        normal_coefs[2],
        interp
    );
}

int mr_fvm_calculate_ghost_cell(mr_ocforest *forest, mr_int coarse_cell_idx, mr_int fine_cell_idx, mr_fvm_interpolation_cb interp) {
    assert(forest);
    assert(coarse_cell_idx != MR_INVALID_INDEX);
    assert(fine_cell_idx != MR_INVALID_INDEX);
    assert(interp.fn);

    mr_octree_cell *coarse_cell = mr_ocforest_get_cell(forest, coarse_cell_idx);
    mr_octree_cell *fine_cell = mr_ocforest_get_cell(forest, fine_cell_idx);
    mr_direction dir = calculate_ghost_cell_direction(coarse_cell, fine_cell);

    mr_int tangent_first_local_idx = get_tangent_stencil_first_local_idx(forest, coarse_cell_idx, coarse_cell, dir);
    mr_int normal_stencil_cell_indices[2] = { 0 };
    get_normal_stencil_indices(forest, fine_cell_idx, dir, normal_stencil_cell_indices);

    mr_float tangent_data_point[MR_NB_AXES] = { 0.0f };
    calculate_tangent_data_point(coarse_cell, fine_cell, dir, tangent_data_point);

    return interpolate_ghost_cell(forest, coarse_cell_idx, tangent_first_local_idx, normal_stencil_cell_indices, tangent_data_point, dir, interp);
}