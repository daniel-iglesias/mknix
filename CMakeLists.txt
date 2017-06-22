project (MkniX)

cmake_minimum_required(VERSION 3.0)
add_definitions(-DHAVE_GMM)
add_definitions(-DHAVE_CUDA)
include(CPack)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -fopenmp --std=c++11 -gdwarf-3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} --std=c++11 -gdwarf-3")

add_subdirectory(src)
add_subdirectory(tests)
