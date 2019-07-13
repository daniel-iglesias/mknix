project (MkniX)

cmake_minimum_required(VERSION 2.8.8)

include(CPack)

if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -w -fopenmp -std=c++11 -gdwarf-3")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fPIC std=c++11 -Wall -Wno-unknown-pragmas -Werror -gdwarf-3")
elseif(WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Wall /DMKNIX_EXPORT /D_SCL_SECURE_NO_DEPRECATE /D_CRT_SECURE_NO_WARNINGS /Qopenmp")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Wall /D_SCL_SECURE_NO_DEPRECATE /D_CRT_SECURE_NO_WARNINGS /DMKNIX_EXPORT")
endif()

add_subdirectory(src)
add_subdirectory(tests)
