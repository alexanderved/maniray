add_library(GLAD STATIC "src/gl.c")
target_include_directories(GLAD PRIVATE "${CMAKE_SOURCE_DIR}/include")