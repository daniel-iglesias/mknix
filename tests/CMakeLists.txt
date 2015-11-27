enable_testing()
#add_executable(testprog main.c)
file( GLOB test1files "${CMAKE_SOURCE_DIR}/tests/input_2triangles.fem.mknix" "2triangles.dat" )
file( COPY ${test1files} DESTINATION ${CMAKE_BINARY_DIR}/tests/ )
set( args "input_2triangles.fem.mknix 2 0" ) 
#message(${CMAKE_COMMAND} ${CMAKE_BINARY_DIR}/src/mknix ${CMAKE_SOURCE_DIR}/tests/${args} )
add_test(runtest1 ${CMAKE_BINARY_DIR}/src/mknix input_2triangles.fem.mknix)
add_test(NAME runtest1diff
    COMMAND ${CMAKE_COMMAND} 
    -Dtest_cmd=${CMAKE_BINARY_DIR}/src/mknix
    -Dtest_args:string=${CMAKE_SOURCE_DIR}/tests/${args}
#    -DSOURCEDIR=${CMAKE_CURRENT_SOURCE_DIR}
    -D output_blessed=${CMAKE_SOURCE_DIR}/tests/2triangfem.mec
    -D output_test=${CMAKE_BINARY_DIR}/tests/2triangfem.mec
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)

