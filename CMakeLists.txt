# This file obviously doesn't work, but we will make it one day, trust
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

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

# Different methods are neede for Linux based systems and Windows
# I don't know anyone with a Mac, so no build available for that

# Unix File ----------------------------------------------

if (UNIX)

    find_package(OpenMP REQUIRED)
    find_package(CUDA QUIET)

    if(CUDA_FOUND)
        set(CMAKE_CUDA_ARCHITECTURES all-major)
        enable_language(CUDA)
        
        set(PROJECT_SOURCES
            src/main.cu
        )
        
        include_directories(
            TPMS-Characterization PUBLIC
            "${PROJECT_SOURCE_DIR}/lib/cuda_solvers"
        )

        set(CMAKE_CUDA_COMPILER "nvcc")

        add_definitions(-DUSE_CUDA)
        
        # Add Executable --------------------------------------------
        add_executable(
            TPMS-Characterization
            ${PROJECT_SOURCES}
        )

        set_property(TARGET TPMS-Characterization PROPERTY CUDA_SEPARABLE_COMPILATION ON)
        set_property(TARGET TPMS-Characterization PROPERTY CUDA_ARCHITECTURES 86)

        target_compile_features(TPMS-Characterization PUBLIC cxx_std_17)

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-fopenmp")
        target_link_libraries(TPMS-Characterization PRIVATE OpenMP::OpenMP_CXX)
    else()
        # CPU version
        set(PROJECT_SOURCES
            src/main.cpp
        )
        # Add Executable --------------------------------------------
        add_executable(
            TPMS-Characterization-CPU
            ${PROJECT_SOURCES}
        )

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        target_link_libraries(TPMS-Characterization-CPU PRIVATE OpenMP::OpenMP_CXX)
        target_compile_features(TPMS-Characterization-CPU PUBLIC cxx_std_17)
    endif(CUDA)

endif (UNIX)

# Windows File (no GPU support) ------------------------------

if(WIN32)
    set(PROJECT_SOURCES
            src/main.cpp
    )
    # Add Executable --------------------------------------------
    add_executable(
        TPMS-Characterization
        ${PROJECT_SOURCES}
    )
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    target_compile_features(TPMS-Characterization PUBLIC cxx_std_17)
endif(WIN32)

# Set Tests Executables -------------------------------------

set(TEST_SOURCES
    tests/saveTPMS.cpp
)

add_executable(
    saveTPMS
    ${TEST_SOURCES}
)

