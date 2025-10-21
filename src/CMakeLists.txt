
set(LMX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LMX-1.0/src)
set(GMM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gmm/gmm-3.0/include)

include_directories(SYSTEM ${LMX_INCLUDE_DIR} ${LMX_INCLUDE_DIR}/LMX)
include_directories(SYSTEM ${GMM_INCLUDE_DIR})

add_subdirectory(core)
add_subdirectory(reader)
add_subdirectory(simulation)
add_subdirectory(system)

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(mknix SHARED ${HEADERS} ${SOURCES}
    $<TARGET_OBJECTS:core>
    $<TARGET_OBJECTS:reader>
    $<TARGET_OBJECTS:simulation>
    $<TARGET_OBJECTS:system>)
    
add_library(mknix_static STATIC ${HEADERS} ${SOURCES}
    $<TARGET_OBJECTS:core>
    $<TARGET_OBJECTS:reader>
    $<TARGET_OBJECTS:simulation>
    $<TARGET_OBJECTS:system>)
