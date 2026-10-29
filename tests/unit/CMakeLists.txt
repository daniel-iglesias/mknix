add_subdirectory(core)
add_subdirectory(simulation)
add_subdirectory(system)

add_executable(mknix_reader_unit_tests
	reader/reader_tests_main.cpp
	reader/reader_test.cpp
	reader/readerconstraints_test.cpp
	reader/readerflex_test.cpp
	reader/readerrigid_test.cpp
	reader/sectionreader_test.cpp)
target_link_libraries(mknix_reader_unit_tests ${VTK_LIBRARIES} mknix)

add_test(NAME mknix_reader_unit_tests COMMAND mknix_reader_unit_tests)