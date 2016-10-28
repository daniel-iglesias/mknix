project (MkniX)

cmake_minimum_required(VERSION 2.8.8)

include(CPack)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -fopenmp --std=gnu++11 -gdwarf-3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -gdwarf-3")

add_subdirectory(src)
add_subdirectory(tests)
