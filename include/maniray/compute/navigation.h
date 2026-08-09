#ifndef _MR_NAVIGATION_H
#define _MR_NAVIGATION_H

#define MR_NB_AXES 3

typedef enum mr_adjacency {
    MR_ADJACENCY_NONE = 0,
    MR_ADJACENCY_FACE = 1,
    MR_ADJACENCY_EDGE = 2,
    MR_ADJACENCY_VERTEX = 3,
} mr_adjacency;

typedef enum mr_axis {
    MR_AXIS_X = 0,
    MR_AXIS_Y = 1,
    MR_AXIS_Z = 2,
} mr_axis;

typedef enum mr_sign {
    MR_SIGN_MINUS = 0,
    MR_SIGN_PLUS = 1,
} mr_sign;

typedef enum mr_direction {
    MR_DIRECTION_MI_X = 0,
    MR_DIRECTION_PL_X = 1,
    MR_DIRECTION_MI_Y = 2,
    MR_DIRECTION_PL_Y = 3,
    MR_DIRECTION_MI_Z = 4,
    MR_DIRECTION_PL_Z = 5,
} mr_direction;

static inline int mr_sign_to_mul(mr_sign sign) {
    return sign * 2 - 1;
}

static inline mr_direction mr_direction_create(mr_axis axis, mr_sign sign) {
    return axis * 2 + sign;
}

static inline mr_axis mr_direction_get_axis(mr_direction dir) {
    return dir / 2;
}

static inline mr_sign mr_direction_get_sign(mr_direction dir) {
    return dir % 2;
}

static inline int mr_direction_get_sign_mul(mr_direction dir) {
    return mr_sign_to_mul(mr_direction_get_sign(dir));
}

static inline mr_direction mr_direction_reflect(mr_direction dir) {
    return dir - mr_direction_get_sign_mul(dir);
}

#endif // _MR_NAVIGATION_H