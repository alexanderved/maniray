add_library(maniray STATIC)

add_subdirectory(display)
add_subdirectory(utils)
add_subdirectory(compute)

target_compile_options(maniray PUBLIC -Wall -Wextra -Wpedantic -Wno-unused-function)
target_compile_options(maniray PUBLIC
    $<$<CONFIG:Debug>:-O0 -g>
    $<$<CONFIG:Release>:-O3 -flto -DNDEBUG>
    $<$<CONFIG:ReleasePerf>:-g -fno-omit-frame-pointer -O2 -DNDEBUG>
)
# target_link_options(maniray PUBLIC -fsanitize=address -fsanitize=leak)

target_link_libraries(maniray ${LIBS})

target_include_directories(maniray PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_include_directories(maniray PUBLIC ${PROJECT_SOURCE_DIR}/external)

set_property(TARGET maniray PROPERTY C_STANDARD 11)