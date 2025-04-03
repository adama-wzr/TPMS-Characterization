# This file obviously doesn't work, but we will make it one day, trust
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# version information--------------------------
project(TPMS-Characterization LANGUAGES CXX CUDA VERSION 1.0)

# Include directories -------------------------
include_directories(
    TPMS-Characterization PUBLIC
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/lib"
    "${PROJECT_SOURCE_DIR}/tests"
    "${PROJECT_SOURCE_DIR}/src_cuda"
)

# Different methods are neede for Linux based systems and Windows
# I don't know anyone with a Mac, so no build available for that

if (WIN32)
    enable_language(CUDA)
endif (WIN32)

if (UNIX)
    find_package(OpenMP REQUIRED)
    find_package(CUDA REQUIRED)
    set(CMAKE_CUDA_COMPILER "nvcc")
endif (UNIX)

# Set Project Sources ----------------------------------------

set(PROJECT_SOURCES
    src/main.cu
    # src/helper.cuh
)

# Add Executable --------------------------------------------

add_executable(
    TPMS-Characterization
    ${PROJECT_SOURCES}
)

# Compilation options to ensure this works with CUDA + OpenMP
set_property(TARGET TPMS-Characterization PROPERTY CUDA_SEPARABLE_COMPILATION ON)
set_property(TARGET TPMS-Characterization PROPERTY CUDA_ARCHITECTURES all-major)
# set_property(TARGET TPMS-Characterization PROPERTY CUDA_SEPARABLE_COMPILATION ON)
target_compile_features(TPMS-Characterization PUBLIC cxx_std_23)

# target_link_libraries(TPMS-Characterization PRIVATE OpenMP::OpenMP_CXX)

if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-fopenmp")

    target_link_libraries(TPMS-Characterization PRIVATE OpenMP::OpenMP_CXX)
endif (UNIX)

# Set Tests Executables -------------------------------------

set(TEST_SOURCES
    tests/saveTPMS.cpp
)

add_executable(
    saveTPMS
    ${TEST_SOURCES}
)

