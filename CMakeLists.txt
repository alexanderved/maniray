cmake_minimum_required(VERSION 3.13)
project(ManiRay)

include(ExternalProject)

if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(OpenMP REQUIRED)

find_package(glfw3 3.3 REQUIRED)
find_package(OpenGL REQUIRED)
set(LIBS m glfw OpenGL::GL GLAD lis OpenMP::OpenMP_C)

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(examples)