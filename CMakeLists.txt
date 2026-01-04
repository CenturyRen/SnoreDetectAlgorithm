cmake_minimum_required(VERSION 3.10)
project(snoreDetectAlgorithm)
set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR}/include)

link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable(snoreDetectAlgorithm src/main.cpp)

target_link_libraries(snoreDetectAlgorithm)