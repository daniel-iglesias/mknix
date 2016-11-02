cmake_minimum_required(VERSION 3.0)

set(CUDA_VERSION 7.5)
find_package(CUDA ${CUDA_VERSION} REQUIRED)

set(CUDA_ARCH -arch=sm_30)
set(CUDA_NVCC_FLAGS ${CUDA_ARCH})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES") #Line needed in ubuntu 16.04. Otherwise complains about memcpy
set(CUDA_INCLUDE_DIRECTORIES /usr/local/cuda-${CUDA_VERSION}/include)
set(CUDA_LIB_DIR /usr/local/cuda-${CUDA_VERSION}/lib64)
set(CUDA_SEPARABLE_COMPILATION ON)
#set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE ON)
set(BUILD_SHARED_LIBS OFF)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CUDA_HOST_COMPILER g++)
set(CUDA_LIBRARIES -L/usr/local/cuda-${CUDA_VERSION}/lib64 cuda cudart)
set(CU nvcc)
#set(DEVICE_LIB device_o)

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)
file(GLOB HEADERS_CUDA *.cuh)
file(GLOB SOURCES_CUDA *.cu)

include_directories(${CUDA_INCLUDE_DIRECTORIES})
include_directories(${PROJECT_SOURCE_DIR}/src)

cuda_add_library(gpu ${SOURCES_CUDA} ${SOURCES})
