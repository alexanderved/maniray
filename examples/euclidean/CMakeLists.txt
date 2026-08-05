add_executable(euclidean euclidean.c)

target_include_directories(euclidean PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(euclidean PRIVATE maniray)
set_property(TARGET euclidean PROPERTY C_STANDARD 11)