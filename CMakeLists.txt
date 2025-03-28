# This file obviously doesn't work, but we will make it one day, trust
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# version information--------------------------
project(TPMS-Characterization VERSION 1.0)

# set(CMAKE_CXX_STANDARD 23)
# set(CMAKE_CXX_STANDARD_REQUIRED True)


# Include directories -------------------------
include_directories("
    ${PROJECT_SOURCE_DIR}

")

enable_language(CUDA)

# find_package(OpenMP)
# find_package(CUDA)



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