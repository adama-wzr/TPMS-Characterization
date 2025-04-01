# TPMS-Characterization
This repository is dedicated to holding code pertaining to the characterization of TPMS structures with finite thicknesses.

# Requiremets

The versions are often soft requirements, meaning they are the only versions we have tested on, but might work on older versions.

- NVIDIA Compute Capability >= 8.6
- CUDA >= 11.5 (tested all the way up to 12.8)
- gcc >= 11.0
- C++17 or newer
- CMAKE 3.0 or newer
- OpenMP (any recent version is fine)

# Compilation

On Ubuntu, assuming CUDA is installed and in the $PATH, the compilation should be simple:

```bash
nvcc -Xcompiler -fopenmp ./main.cu
```
Will produce your executable `a.out`. Electively use the flag `-o anyName` to give a name for the executable.

On Windows, it some more flags might be necessary to control execution. The main issues are with different C++ versions (or not having a default) and sometimes the NVIDIA compiler has trouble identifying the architecture of the GPU. To compile on Windows, use the following command:

```bash
nvcc -Xcompiler -openmp -std=c++17 -arch=sm_XX .\main.cu
```

In the `-arch=sm_XX`, replace the `XX` with your GPU compute capability (check it (here)[https://developer.nvidia.com/cuda-gpus]).

**New Method**: Instructions for how to build the code:

- Download all the files from this repository. 
- Create a new folder named "build".
- Enter the new folder.
- Use the following command: `cmake ..`
- If all requirements are met, then the command above should be successful
- Run `cmake --build . --config Release`
- This will create a folder "Release" inside of the "build" folder.
- Inside the "Release" folder, the executable "TPMS-Executable.exe" can be found (Windows).

The current CMAKE file does not work on Ubuntu/Rocky (HPC), it fails to link the OpenMP library with the CUDA calls. I have no idea why it works on Windows but not on Linux.

# Acknowledgements

This work wouldn't be possible without the computational time awarded as part of the following grants:

This work used Expanse(GPU) at SDSC through allocations MAT210014 and MAT230071 from the Advanced Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS) program, which is supported by National Science Foundation grants #2138259, #2138286, #2138307, #2137603, and #2138296.

# Upcoming Changes
