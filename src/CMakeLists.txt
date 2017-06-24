cmake_minimum_required(VERSION 3.0)

add_definitions(-DHAVE_GMM)
add_definitions(-DHAVE_CUDA)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -fopenmp --std=c++11 -gdwarf-3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --std=c++11 -ggdb")
set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)
set(GPU_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src/gpu)

#if(HAVE_CUDA)
set(CUDA_INCLUDE_DIRECTORIES /usr/local/cuda-8.0/include)
#endif()

include_directories(${LMX_INCLUDE_DIR})
include_directories(${GMM_INCLUDE_DIR})
include_directories(${GPU_INCLUDE_DIR})
#if(HAVE_CUDA)
include_directories(${CUDA_INCLUDE_DIRECTORIES})
#endif()

add_subdirectory(core)
add_subdirectory(reader)
add_subdirectory(simulation)
add_subdirectory(system)
add_subdirectory(cpu)
add_subdirectory(gpu)

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

#if(HAVE_CUDA)
add_library(mknix ${HEADERS} ${SOURCES}
    $<TARGET_OBJECTS:core>
    $<TARGET_OBJECTS:reader>
    $<TARGET_OBJECTS:simulation>
    $<TARGET_OBJECTS:system>
    $<TARGET_OBJECTS:cpu>)

#link_directories(${PROJECT_SOURCE_DIR}/src/gpu)
target_link_libraries(mknix gpu)
#else()
#add_library(mknix ${HEADERS} ${SOURCES}
#    $<TARGET_OBJECTS:core>
#    $<TARGET_OBJECTS:reader>
#    $<TARGET_OBJECTS:simulation>
#    $<TARGET_OBJECTS:system>
#    $<TARGET_OBJECTS:gpu>)
#endif()
