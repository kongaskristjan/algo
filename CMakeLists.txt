cmake_minimum_required(VERSION 3.10)
project(algo)

set(CMAKE_CXX_STANDARD 14)

find_package(GTest REQUIRED)
add_executable(RunTests Tests.cpp algo.h Tests.cpp)
target_link_libraries(RunTests GTest::GTest GTest::Main)

gtest_discover_tests(RunTests)
add_test(NAME monolithic COMMAND RunTests)
