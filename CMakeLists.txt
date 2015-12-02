project (MkniX)

cmake_minimum_required(VERSION 2.8.8)

include(CPack)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -Werror -fopenmp -std=gnu++0x -gdwarf-3")

add_subdirectory(src)
add_subdirectory(tests)


