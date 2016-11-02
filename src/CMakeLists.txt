cmake_minimum_required(VERSION 3.0)

#find_package(CUDA REQUIRED)
#set(CUDA_ARCH -arch=sm_30)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${CUDA_ARCH} ; -Xcompiler; -fPIC; --std=c++11; -rdc=true)
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${CUDA_ARCH} ; -Xcompiler; -fPIC; --std=c++11; -rdc=true)
set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)
set(GPU_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src/gpu)
set(CUDA_INCLUDE_DIRECTORIES /usr/local/cuda-7.5/include)
#set(CUDA_SEPARABLE_COMPILATION ON)
#set(CUDA_PROPAGATE_HOST_FLAGS OFF)
#set(CUDA_HOST_COMPILER clang++)
#set(CUDA_LIBRARIES cudart)
#set(CU nvcc)
#set(DEVICE_LIB device_o)

add_definitions(-DHAVE_GMM)
include_directories(${LMX_INCLUDE_DIR})
include_directories(${GMM_INCLUDE_DIR})
include_directories(${GPU_INCLUDE_DIR})
include_directories(${CUDA_INCLUDE_DIRECTORIES})

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
    $<TARGET_OBJECTS:system>)

#link_directories( ${PROJECT_SOURCE_DIR}/src/gpu)
target_link_libraries(mknix gpu)
