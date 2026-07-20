cmake_minimum_required(VERSION 3.13)
project(maniray)

if(POLICY CMP0072)
    cmake_policy(SET CMP0072 NEW)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(glfw3 3.3 REQUIRED)
find_package(OpenGL REQUIRED)
set(LIBS m glfw OpenGL::GL GLAD)

add_subdirectory(external)
add_subdirectory(src)

target_include_directories(maniray PRIVATE include)
set_property(TARGET maniray PROPERTY C_STANDARD 11)