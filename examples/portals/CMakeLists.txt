add_executable(portals portals.c)

target_include_directories(portals PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(portals PRIVATE maniray)
set_property(TARGET portals PROPERTY C_STANDARD 11)