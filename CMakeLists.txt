cmake_minimum_required(VERSION 3.12)

project(fastwindingnumber LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(fastwindingnumber
    src/FastWindingNumber.cpp
    src/IO.cpp
)

target_include_directories(fastwindingnumber
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if (MSVC)
    target_compile_options(fastwindingnumber PRIVATE /W4)
else()
    target_compile_options(fastwindingnumber PRIVATE -Wall -Wextra -pedantic)
endif()

add_executable(fwn_query tools/fwn_query.cpp)
target_link_libraries(fwn_query PRIVATE fastwindingnumber)

add_executable(fwn_benchmark tools/fwn_benchmark.cpp)
target_link_libraries(fwn_benchmark PRIVATE fastwindingnumber)

add_executable(fwn_basic_test tests/basic_test.cpp)
target_link_libraries(fwn_basic_test PRIVATE fastwindingnumber)

enable_testing()
add_test(NAME fwn_basic_test COMMAND fwn_basic_test)
