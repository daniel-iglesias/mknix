cmake_minimum_required(VERSION 3.0)

enable_testing()

add_definitions(-DHAVE_GMM)
set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${LMX_INCLUDE_DIR})
include_directories(${GMM_INCLUDE_DIR})

add_executable(mknixrunner mknix.cpp)

#link_directories( ${PROJECT_SOURCE_DIR}/src/gpu)
target_link_libraries(mknixrunner mknix)

set(test1dir ${CMAKE_CURRENT_BINARY_DIR}/test1)

set(test1files input_2triangles.fem.mknix 2triangles.dat)
file(COPY ${test1files} DESTINATION ${test1dir})

set(args input_2triangles.fem.mknix 2 0)

# need to replace ; with \\; so that args gets passed as semi colon separated list to command
string(REPLACE ";" "\\;" args "${args}")

add_test(NAME test1
    COMMAND ${CMAKE_COMMAND}
      -Dtest_cmd=${CMAKE_CURRENT_BINARY_DIR}/mknixrunner
      -Dtest_args=${args}
      -Doutput_blessed=${CMAKE_SOURCE_DIR}/tests/2triangfem.mec
      -Doutput_test=2triangfem.mec
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    WORKING_DIRECTORY ${test1dir})
