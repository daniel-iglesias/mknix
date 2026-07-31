enable_testing()

set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${LMX_INCLUDE_DIR})
include_directories(${GMM_INCLUDE_DIR})

add_executable(mknixrunner mknix.cpp)
target_link_libraries(mknixrunner ${VTK_LIBRARIES} mknix)

add_executable(mknixwrapper mknixwrapper.cpp)
target_link_libraries(mknixwrapper ${VTK_LIBRARIES} mknix )

include(mknix_test.cmake)

MKNIX_TEST(
    NAME test1_thermalfluence
    COMMAND mknixrunner
    COMMAND_ARGS input_2triangles.fem.mknix 2 0
    INPUT_FILES input_2triangles.fem.mknix 2triangles.dat
    OUTPUT_FILE 2triangfem.mec
)

MKNIX_TEST(
    NAME test2_thermalfluence_timefile
    COMMAND mknixrunner
    COMMAND_ARGS input_2triangles_timefile.fem.mknix 2 0
    INPUT_FILES input_2triangles_timefile.fem.mknix 2triangles.dat timescale.dat
    OUTPUT_FILE 2triangfem_timefile.mec
)
