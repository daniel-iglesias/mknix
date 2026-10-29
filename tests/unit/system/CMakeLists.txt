add_executable(mknix_system_unit_tests
    system_tests_main.cpp
    system_motion_test.cpp
    load_constraint_test.cpp)
target_link_libraries(mknix_system_unit_tests ${VTK_LIBRARIES} mknix)

add_test(NAME mknix_system_unit_tests COMMAND mknix_system_unit_tests)