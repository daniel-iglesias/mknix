cmake_minimum_required(VERSION 2.8.8)

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

include_directories(${PROJECT_SOURCE_DIR}/src)

add_library(core OBJECT ${SOURCES} ${HEADERS})