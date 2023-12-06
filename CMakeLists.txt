cmake_minimum_required(VERSION 3.27)

project(MainScf)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 17)

add_executable(MainScf main_scf.cpp)

# add_subdirectory follows the CMakeLists.txt in the specified directory 
add_subdirectory(src/)
add_subdirectory(tests/)
add_subdirectory(third_party/googletest/)

# src is the name of library (aka. consistent with src in target_include_directories(src PUBLIC <include directories>))
target_link_libraries(MainScf PUBLIC src)

enable_testing ()
# NAME: name of the test 
# COMMAND: existing executable 
add_test (NAME ScfTest COMMAND ScfTest)