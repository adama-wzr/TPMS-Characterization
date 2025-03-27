# TPMS-Characterization
This repository is dedicated to holding code pertaining to the characterization of TPMS structures with finite thicknesses.

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

# Requiremets

# Acknowledgements

# Upcoming Changes
