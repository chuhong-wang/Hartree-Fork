cmake_minimum_required(VERSION 3.27)

set(This Molecule_1)

project(${This} C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

enable_testing()

# include 
include_directories(include)
include_directories(third_party/)

# source file 
file(GLOB SOURCES "src/*.cpp")

add_library(${This} STATIC ${SOURCES})

# Third-party libraries (googletest and Eigen)
add_subdirectory(third_party/googletest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# Test execuatable
file(GLOB TEST_SOURCES "tests/*.cpp") 
add_executable(MoleculeTest ${TEST_SOURCES})

# Link Google Test and source code to the test executable 
target_link_libraries(MoleculeTest gtest gtest_main ${This})

file(COPY data DESTINATION ${CMAKE_BINARY_DIR})

add_test(
    NAME MoleculeTest
    COMMAND MoleculeTest
)