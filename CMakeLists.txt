# This file obviously doesn't work, but we will make it one day, trust
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)



# set(CMAKE_CXX_STANDARD 23)
# set(CMAKE_CXX_COMPILER "gcc")
set(CMAKE_CUDA_COMPILER "nvcc")

# version information--------------------------
project(TPMS-Characterization LANGUAGES CXX CUDA VERSION 1.0)

# Include directories -------------------------
include_directories("
    ${PROJECT}
")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-fopenmp")

# enable_language(CUDA)

find_package(OpenMP REQUIRED)
find_package(CUDA REQUIRED)

set(PROJECT_SOURCES
    src/main.cu
    src/helper.cuh
)

add_executable(
    TPMS-Characterization
    ${PROJECT_SOURCES}
)

set_property(TARGET TPMS-Characterization PROPERTY CUDA_ARCHITECTURES all-major)
target_compile_features(TPMS-Characterization PUBLIC cxx_std_23)
target_link_libraries(TPMS-Characterization PRIVATE OpenMP::OpenMP_CXX)