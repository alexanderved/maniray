target_sources(maniray PRIVATE manifold.c
    octree.c
    math.c
    codes.c
    manifold3d.c)

add_subdirectory(fvm)
add_subdirectory(linear_system)