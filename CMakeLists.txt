cmake_minimum_required(VERSION 3.27)

project(MainScf)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED 17)

# include_directories(third_party/)
# include_directories(include/)

add_executable(MainScf main_scf.cpp)


add_subdirectory(src/)
add_subdirectory(tests/)
add_subdirectory(third_party/googletest/)

target_link_libraries(MainScf PUBLIC src)

enable_testing ()
add_test (NAME ScfTest COMMAND ScfTest)