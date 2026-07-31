enable_testing()

set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${LMX_INCLUDE_DIR})
include_directories(${GMM_INCLUDE_DIR})

add_executable(mknixrunner mknix.cpp)
target_link_libraries(mknixrunner ${VTK_LIBRARIES} mknix)

add_executable(mknixwrapper mknixwrapper.cpp)
target_link_libraries(mknixwrapper ${VTK_LIBRARIES} mknix)

set(MKNIX_RUNNER_PATH     ${CMAKE_CURRENT_BINARY_DIR}/mknixrunner)
set(MKNIX_RESOURCES_DIR   ${CMAKE_CURRENT_SOURCE_DIR}/resources)
set(MKNIX_RUN_TEST_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake)

include(mknix_test.cmake)

# Test cases in alphabetical order by directory name
add_subdirectory(test1_thermalfluence)
add_subdirectory(test2_thermalfluence_timefile)
add_subdirectory(test3_thermalflux_timefile)
