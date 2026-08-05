add_library(maniray STATIC)

add_subdirectory(display)
add_subdirectory(utils)
add_subdirectory(compute)

target_compile_options(maniray PUBLIC -Wall -Wextra -Wpedantic)
# target_link_options(maniray PRIVATE -fsanitize=address -fsanitize=leak)

target_link_libraries(maniray ${LIBS})

target_include_directories(maniray PUBLIC ${PROJECT_SOURCE_DIR}/include)
set_property(TARGET maniray PROPERTY C_STANDARD 11)