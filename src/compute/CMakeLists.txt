target_sources(maniray PRIVATE manifold.c
    octree.c
    sparse_matrix.c
    math.c
    codes.c
    manifold3d.c)

add_subdirectory(fvm)