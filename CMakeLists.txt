# This file obviously doesn't work, but we will make it one day, trust
cmake_minimum_required(VERSION 3.00)

find_package(OpenMP)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(TPMS_Characterization VERSION 1.0)

add_executable(
    TPMS_Characterization main.cu helper.cuh
)
