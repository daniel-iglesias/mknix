add_executable(mknix_simulation_unit_tests
    simulation_tests_main.cpp
    analysis_test.cpp
    simulation_settings_test.cpp)
target_link_libraries(mknix_simulation_unit_tests ${VTK_LIBRARIES} mknix)

add_test(NAME mknix_simulation_unit_tests COMMAND mknix_simulation_unit_tests)