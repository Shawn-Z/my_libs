cmake_minimum_required(VERSION 3.12)
project(my_libs)

set(CMAKE_CXX_STANDARD 11)

include_directories(
    shawn
)

add_executable(my_libs
    shawn/test/main.cpp
    shawn/STime.hpp
    shawn/SVector.hpp
    shawn/SJerk.hpp
    shawn/SPoints.hpp
    shawn/SThreePointsCurvature.hpp
)
