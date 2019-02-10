cmake_minimum_required(VERSION 3.0.0)
project(Graphs VERSION 0.1.0)

include(CTest)
enable_testing()

add_executable(Graphs main.cpp Graph.cpp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
