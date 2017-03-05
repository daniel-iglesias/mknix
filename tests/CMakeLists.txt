cmake_minimum_required(VERSION 2.8.8)

enable_testing()

set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${LMX_INCLUDE_DIR})
include_directories(${GMM_INCLUDE_DIR})

add_executable(mknixrunner mknix.cpp)
target_link_libraries(mknixrunner mknix)

include(mknix_test.cmake)

MKNIX_TEST(
    NAME test1
    COMMAND_ARGS input_2triangles.fem.mknix 2 0
    INPUT_FILES input_2triangles.fem.mknix 2triangles.dat
    OUTPUT_FILE 2triangfem.mec
)
