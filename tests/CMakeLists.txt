cmake_minimum_required(VERSION 3.10)

include(GoogleTest)

# include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
file(GLOB TEST_SOURCES "${PROJECT_SOURCE_DIR}/tests/*.cpp")
# Create the test executable
add_executable(ScfTest ${TEST_SOURCES})


# Link Google Test and source code to the test executable 
target_link_libraries(ScfTest gtest gtest_main src)

# copy data/ to tests/ so that the tests are executed in the same architecture as the main programme 
file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR}/tests)


