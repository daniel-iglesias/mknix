cmake_minimum_required(VERSION 2.8.8)

set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)

add_definitions(-DHAVE_GMM)
include_directories(${LMX_INCLUDE_DIR})
include_directories(${GMM_INCLUDE_DIR})

add_subdirectory(core)
add_subdirectory(reader)
add_subdirectory(simulation)
add_subdirectory(system)

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(mknix ${HEADERS} ${SOURCES}
    $<TARGET_OBJECTS:core>
    $<TARGET_OBJECTS:reader>
    $<TARGET_OBJECTS:simulation>
    $<TARGET_OBJECTS:system>)