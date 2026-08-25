add_subdirectory(glad)

set(LIS_INSTALL_DIR "${CMAKE_BINARY_DIR}/external/lis")
ExternalProject_Add(
    lis_external
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/lis"
    PREFIX            "${LIS_INSTALL_DIR}"
    CONFIGURE_COMMAND  <SOURCE_DIR>/configure --prefix=${LIS_INSTALL_DIR} --enable-omp
    BUILD_COMMAND      make
    INSTALL_COMMAND    make install
)

add_library(lis STATIC IMPORTED GLOBAL)
set_target_properties(
    lis PROPERTIES
    IMPORTED_LOCATION ${LIS_INSTALL_DIR}/lib/liblis.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/lis/include
)
add_dependencies(lis lis_external)