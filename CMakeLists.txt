cmake_minimum_required(VERSION 3.14)
project(reax_tools)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

file(GLOB_RECURSE SRC_FILES src/cpp/*.cpp)
list(APPEND SRC_FILES
    "${CMAKE_SOURCE_DIR}/src/cpp/fmt/format.cc"
    "${CMAKE_SOURCE_DIR}/src/cpp/fmt/os.cc"
)

add_executable(reax_tools_core ${SRC_FILES})

target_include_directories(reax_tools_core PRIVATE
    src/cpp
    src/cpp/fmt
)
