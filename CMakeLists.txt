cmake_minimum_required(VERSION 3.9)
project(integral)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Og")

add_executable(integral main.cpp)
target_link_libraries(integral ${CMAKE_THREAD_LIBS_INIT})