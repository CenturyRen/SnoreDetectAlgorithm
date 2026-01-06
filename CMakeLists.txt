cmake_minimum_required(VERSION 3.10)
project(snoreDetectAlgorithm)
set(CMAKE_CXX_STANDARD 17)

include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/m)
link_directories(${CMAKE_SOURCE_DIR}/lib)

add_executable(snoreDetectAlgorithm 
                src/main.cpp
                src/buildMFCC.cpp
                src/detectSnore.cpp)

target_link_libraries(snoreDetectAlgorithm)