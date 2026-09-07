#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "maniray/compute/math.h"
#include "maniray/compute/fvm/interpolation.h"
#include "maniray/compute/fvm/cell.h"
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

static mr_direction calculate_fine_cell_normal_direction(mr_octree_cell *coarse_cell, mr_octree_cell *fine_cell) {
    for (mr_axis axis = MR_AXIS_X; axis <= MR_AXIS_Z; ++axis) {
        mr_float diff = get_cell_center_axis(fine_cell, axis) - get_cell_center_axis(coarse_cell, axis);
        if (MR_ABS(diff) > coarse_cell->dim / 2.0f) {
            mr_sign sign = diff < 0.0 ? MR_SIGN_MINUS : MR_SIGN_PLUS;
            // In case of periodic boundaries cells might be located at opposite sides;
            // therefore, the sign must be changed to the opposite one when the distance
            // between cells is larger than the dimension of the coarse cell.
            mr_sign periodic_sign = MR_ABS(diff) < coarse_cell->dim ? sign : (1 - sign);

            return mr_direction_create(axis, periodic_sign);
        }
    }

    assert(false);
    return MR_DIRECTION_MI_X;
}

static void get_normal_stencil_indices(mr_ocforest *forest, mr_int fine_cell_idx, mr_direction dir, mr_int stencil[2]) {
    mr_octree_cell_neighbor n = mr_octree_find_face_neighbor_cells(forest, fine_cell_idx, dir);
    assert(n.type == MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE);

    stencil[0] = n.neighbor_idx;
    stencil[1] = fine_cell_idx;
}

static void calculate_normal_coefs(mr_float interp_coord, mr_float coefs[Q_STENCIL_DIM]) {
    mr_float stencil[Q_STENCIL_DIM] = { 0.0f, 1.0f, 2.5f };

    for (size_t i = 0; i < Q_STENCIL_DIM; ++i) {
        coefs[i] = lagrange_coef(interp_coord, stencil, i);
    }
}

static int interpolate_coarse_fine_normal(
    mr_ocforest *forest,
    mr_int fine_cell_idx,
    mr_direction dir,
    mr_float interp_coord,
    mr_fvm_interpolation_cb interp,
    mr_float *tangent_mul
) {
    mr_int normal_stencil_cell_indices[2] = { 0 };
    get_normal_stencil_indices(forest, fine_cell_idx, dir, normal_stencil_cell_indices);

    mr_float normal_coefs[Q_STENCIL_DIM] = { 0.0f };
    calculate_normal_coefs(interp_coord, normal_coefs);

    for (size_t i = 0; i < 2; ++i) {
        mr_int cell_idx = normal_stencil_cell_indices[i];
        mr_float coef = normal_coefs[i];

        if (interp.fn(forest, cell_idx, coef, interp.userdata) != MR_SUCCESS) {
            return MR_FAILURE;
        }
    }

    *tangent_mul = normal_coefs[2];

    return MR_SUCCESS;
}

static void get_tangent_axes(mr_axis normal_axis, mr_axis tangent_axes[2]) {
    tangent_axes[0] = normal_axis == MR_AXIS_X ? MR_AXIS_Y : MR_AXIS_X;
    tangent_axes[1] = normal_axis == MR_AXIS_Z ? MR_AXIS_Y : MR_AXIS_Z;
}

static void get_host_point_tangent_coords(
    mr_ocforest *forest,
    mr_int coarse_cell_idx,
    mr_octree_cell *coarse_cell,
    mr_axis normal_axis,
    mr_int coords[2]
) {
    mr_octree_node *node = mr_ocforest_get_node(forest, coarse_cell->parent);
    mr_int local_idx = coarse_cell_idx - node->first_child;

    mr_int i = 0;
    for (mr_axis axis = MR_AXIS_X; axis <= MR_AXIS_Z; ++axis) {
        if (axis == normal_axis) {
            continue;
        }

        mr_int local_coord = extract_coord(local_idx, axis);
        coords[i++] = (local_coord - get_first_cell_local_coord(local_coord)) * 4;
    }
}

static void get_tangent_interp_point_offset(
    mr_octree_cell *coarse_cell,
    mr_octree_cell *fine_cell,
    mr_axis normal_axis,
    mr_int offset[2]
) {
    mr_int i = 0;
    for (mr_axis axis = MR_AXIS_X; axis <= MR_AXIS_Z; ++axis) {
        if (axis == normal_axis) {
            continue;
        }

        mr_float diff = get_cell_center_axis(fine_cell, axis) - get_cell_center_axis(coarse_cell, axis);
        offset[i++] = diff > 0.0 ? 1 : -1;
    }
}

static mr_int get_tangent_stencil_first_local_idx(
    mr_ocforest *forest,
    mr_int coarse_cell_idx,
    mr_octree_cell *cell,
    mr_axis normal_axis
) {
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);
    mr_int local_idx = coarse_cell_idx - node->first_child;

    mr_int local_coords[] = {
        extract_coord(local_idx, MR_AXIS_X),
        extract_coord(local_idx, MR_AXIS_Y),
        extract_coord(local_idx, MR_AXIS_Z),
    };

    mr_int first_local_coords[MR_NB_AXES] = { 0 };
    for (mr_axis axis = MR_AXIS_X; axis <= MR_AXIS_Z; ++axis) {
        first_local_coords[axis] = (axis == normal_axis) ? local_coords[axis]
                                                         : get_first_cell_local_coord(local_coords[axis]);
    }

    return pack_local_idx(first_local_coords);
}

static int interpolate_coarse_fine_tangent(
    mr_ocforest *forest,
    mr_int coarse_cell_idx,
    mr_int fine_cell_idx,
    mr_direction dir,
    mr_float mul,
    mr_fvm_interpolation_cb interp
) { 
    mr_octree_cell *coarse_cell = mr_ocforest_get_cell(forest, coarse_cell_idx);
    mr_octree_cell *fine_cell = mr_ocforest_get_cell(forest, fine_cell_idx);

    mr_axis normal_axis = mr_direction_get_axis(dir);
    mr_int tangent_first_local_idx = get_tangent_stencil_first_local_idx(forest, coarse_cell_idx, coarse_cell, normal_axis);

    mr_int host_point_coords[2] = { 0 };
    get_host_point_tangent_coords(forest, coarse_cell_idx, coarse_cell, normal_axis, host_point_coords);

    mr_int interp_offset[2] = { 0 };
    get_tangent_interp_point_offset(coarse_cell, fine_cell, normal_axis, interp_offset);

    mr_int interp_coords[] = { host_point_coords[0] + interp_offset[0], host_point_coords[1] + interp_offset[1] };
    mr_octree_node *node = mr_ocforest_get_node(forest, coarse_cell->parent);

    mr_axis tangent_axes[2] = { 0 };
    get_tangent_axes(normal_axis, tangent_axes);

    for (mr_int i = 0; i < MR_NB_AXES; ++i) {
        for (mr_int j = 0; j < MR_NB_AXES; ++j) {
            mr_int local_idx = tangent_first_local_idx + (i << tangent_axes[0] * 2) + (j << tangent_axes[1] * 2);
            mr_int cell_idx = node->first_child + local_idx;

            mr_float coef = mul
                * lagrange_coef_with_idx(interp_coords[0] / 4.0f, 0, i)
                * lagrange_coef_with_idx(interp_coords[1] / 4.0f, 0, j);
            if (interp.fn(forest, cell_idx, coef, interp.userdata) != MR_SUCCESS) {
                return MR_FAILURE;
            }
        }
    }

    return MR_SUCCESS;
}

static int interpolate_coarse_fine(
    mr_ocforest *forest,
    mr_int coarse_cell_idx,
    mr_int fine_cell_idx,
    mr_float interp_normal_coord,
    mr_fvm_interpolation_cb interp
) {
    mr_octree_cell *coarse_cell = mr_ocforest_get_cell(forest, coarse_cell_idx);
    mr_octree_cell *fine_cell = mr_ocforest_get_cell(forest, fine_cell_idx);
    mr_direction dir = calculate_fine_cell_normal_direction(coarse_cell, fine_cell);

    mr_float tangent_mul = 0.0f;
    int res = interpolate_coarse_fine_normal(forest, fine_cell_idx, dir, interp_normal_coord, interp, &tangent_mul);
    if (res != MR_SUCCESS) {
        return res;
    }

    return interpolate_coarse_fine_tangent(
        forest,
        coarse_cell_idx,
        fine_cell_idx,
        dir,
        tangent_mul,
        interp
    );
}

int mr_fvm_calculate_ghost_cell(mr_ocforest *forest, mr_int coarse_cell_idx, mr_int fine_cell_idx, mr_fvm_interpolation_cb interp) {
    assert(forest);
    assert(coarse_cell_idx != MR_INVALID_INDEX);
    assert(fine_cell_idx != MR_INVALID_INDEX);
    assert(interp.fn);

    return interpolate_coarse_fine(forest, coarse_cell_idx, fine_cell_idx, 2.0f, interp);
}

int mr_fvm_interpolate_face_value(mr_ocforest *forest, mr_int cell_idx, mr_direction face_dir, mr_fvm_interpolation_cb interp) {
    assert(forest);
    assert(cell_idx != MR_INVALID_INDEX);
    assert(interp.fn);

    mr_octree_cell *cell = mr_ocforest_get_cell(forest, cell_idx);
    mr_octree_node *node = mr_ocforest_get_node(forest, cell->parent);

    mr_octree_cell_neighbor neighbor = mr_octree_find_face_neighbor_cells(forest, cell_idx, face_dir);

    switch (neighbor.type) {
        case MR_OCTREE_CELL_NEIGHBOR_NONE: ;
            mr_float face_center[MR_NB_AXES] = { 0.0f };
            mr_cell_face_center(forest, cell_idx, face_dir, face_center);

            return mr_fvm_perform_interpolation(forest, node->root, face_center, interp);

        case MR_OCTREE_CELL_NEIGHBOR_EQUAL_SIZE:
            if (interp.fn(forest, cell_idx, 0.5, interp.userdata) != MR_SUCCESS) {
                return MR_FAILURE;
            }
            return interp.fn(forest, neighbor.neighbor_idx, 0.5, interp.userdata);
        
        default:
            // TODO: Implement for coarse-fine
            abort();
    }
}