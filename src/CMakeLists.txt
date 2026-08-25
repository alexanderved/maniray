add_library(maniray STATIC
    maniray.c)

add_subdirectory(display)
add_subdirectory(utils)
add_subdirectory(compute)

target_compile_options(maniray PUBLIC -Wall -Wextra -Wpedantic -Wno-unused-function)
target_compile_options(maniray PUBLIC
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3 -flto -march=native -DNDEBUG>
    $<$<CONFIG:ReleasePerf>:-g -fno-omit-frame-pointer -O2 -DNDEBUG>
)
# target_link_options(maniray PUBLIC -fsanitize=address -fsanitize=leak)

target_link_libraries(maniray PUBLIC ${LIBS})

target_include_directories(maniray PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(maniray PUBLIC ${PROJECT_SOURCE_DIR}/external)
target_include_directories(maniray PUBLIC ${PROJECT_SOURCE_DIR}/external/lis/include)

set_property(TARGET maniray PROPERTY C_STANDARD 11)