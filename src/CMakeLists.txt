add_executable(maniray main.c)
add_subdirectory(display)
add_subdirectory(utils)
add_subdirectory(compute)

target_compile_options(maniray PRIVATE -Wall -Wextra -Wpedantic)
# target_link_options(maniray PRIVATE -fsanitize=address -fsanitize=leak)

target_link_libraries(maniray ${LIBS})