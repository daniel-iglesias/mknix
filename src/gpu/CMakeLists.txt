cmake_minimum_required(VERSION 2.8.8)

find_package(CUDA REQUIRED)
file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)
file(GLOB HEADERS_CUDA *.cuh)
file(GLOB SOURCES_CUDA *.cu)
CUDA_COMPILE(CU_O ${SOURCES_CUDA} ${HEADERS_CUDA})

include_directories(${PROJECT_SOURCE_DIR}/src)

#add_library(gpu_host OBJECT ${SOURCES} ${HEADERS})
#add_library(gpu_device OBJECT ${SOURCES_CUDA} ${HEADERS_CUDA})
add_library(cuda_compile OBJECT ${SOURCES} ${HEADERS} ${SOURCES_CUDA} ${HEADERS_CUDA})
