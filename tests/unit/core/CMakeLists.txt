add_executable(mknix_core_unit_tests
    core_tests_main.cpp
    material_test.cpp
    point_node_test.cpp
    shape_function_test.cpp
    cell_gausspoint_test.cpp)
target_link_libraries(mknix_core_unit_tests ${VTK_LIBRARIES} mknix)

add_test(NAME mknix_core_unit_tests COMMAND mknix_core_unit_tests)