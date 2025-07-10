# This file obviously doesn't work, but we will make it one day, trust
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# version information--------------------------
project(TPMS-Characterization LANGUAGES CXX VERSION 1.0)

# Include directories -------------------------
include_directories(
    TPMS-Characterization PUBLIC
    "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/lib"
    "${PROJECT_SOURCE_DIR}/tests"
    "${PROJECT_SOURCE_DIR}/lib/cpu_solvers"
)

if (UNIX)
    find_package(OpenMP REQUIRED)
endif (UNIX)

# Different methods are neede for Linux based systems and Windows
# I don't know anyone with a Mac, so no build available for that

# Look for CUDA ------------------------------------------

find_package(CUDA QUIET)

if(CUDA_FOUND)

    include_directories(
        TPMS-Characterization PUBLIC
        "${PROJECT_SOURCE_DIR}/lib/cuda_solvers"
    )

    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE CUDA_ARCHITECTURES all-major)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-fopenmp")

    add_definitions(-DUSE_CUDA)

    if (UNIX)
        set(CMAKE_CUDA_COMPILER "nvcc")
    endif()

    set(PROJECT_SOURCES
        src/main.cu
    )

else()

    set(PROJECT_SOURCES
        src/main.cpp
    )

endif()

# Set Project Sources ----------------------------------------

# set(PROJECT_SOURCES
#     src/main.cpp
# )

# Add Executable --------------------------------------------

add_executable(
    TPMS-Characterization
    ${PROJECT_SOURCES}
)

# Compilation options to ensure this works with CUDA + OpenMP
# if(CUDA_FOUND)
#     set_property(TARGET TPMS-Characterization PROPERTY CUDA_SEPARABLE_COMPILATION ON)
#     set_property(TARGET TPMS-Characterization PROPERTY CUDA_ARCHITECTURES all-major)
#     set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-fopenmp")
# endif()

target_compile_features(TPMS-Characterization PUBLIC cxx_std_17)

if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
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

