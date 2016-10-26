cmake_minimum_required(VERSION 2.8.8)

find_package(CUDA REQUIRED)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_30,code=sm_30; -std=c++11)
set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)
set(GPU_INCLUDE_DIR /usr/local/cuda-7.5/include)
set(CUDA_SEPARABLE_COMPILATION ON)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_HOST_COMPILER clang++)

add_definitions(-DHAVE_GMM)
include_directories(${LMX_INCLUDE_DIR})
include_directories(${GMM_INCLUDE_DIR})
include_directories(${GPU_INCLUDE_DIR})

add_subdirectory(core)
add_subdirectory(reader)
add_subdirectory(simulation)
add_subdirectory(system)
add_subdirectory(gpu)

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)


add_library(mknix ${HEADERS} ${SOURCES}
    $<TARGET_OBJECTS:core>
    $<TARGET_OBJECTS:reader>
    $<TARGET_OBJECTS:simulation>
    $<TARGET_OBJECTS:system>
    $<TARGET_OBJECTS:gpu_host>
    $<TARGET_OBJECTS:gpu_device>)
