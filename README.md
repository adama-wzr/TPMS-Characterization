# TPMS-Characterization

This package contains a collection of code for generating voxelated Triply Periodic Minimal Surfaces (TPMS) and simulating their effective properties. The code is built in a modular and flexible way, since the main goal of this project is to study and classify TPMS structures.

# Table of Contents

1. [Requirements](#requirements)
2. [Compilation](#compilation)
3. [Publications](#publications)
4. [Authors](#authors)
5. [Acknowledgements](#acknowledgements)
6. [Upcoming Changes](#upcoming-changes)
7. [References](#references)

# Requiremets

The versions are often soft requirements, meaning they are the only versions we have tested on, but might work on older versions.

- gcc >= 11.0
- C++17 or newer
- CMAKE 3.15 or newer
- OpenMP (any recent version is fine)

**Soft Requirement**

- NVIDIA GPU (only tested with compute capability >= 7.0)
- CUDA >= 11.5 (tested all the way up to 12.8)

The tortuosity calculation on GPU is hundreds of times faster than CPU and supports multi-GPU execution. GPU is highly recommended for evaluating properties of large (> 300 voxels on each side) structures.

# Compilation

On Windows, to avoid issues, just have cmake and all other requirements installed. Use Visual Studio (the IDE, not VS Code) to build the project for you. In some cases, permission issues might arise when compiling directly on the PowerShell.

Instructions for how to build the code:

- Download all the files from this repository. 
- Create a new folder named "build".
- Enter the new folder.
- Run cmake from inside the build folder: `cmake ..`
- If all requirements are met, then the command above should be successful
- Run `cmake --build . --config Release`
- Windows:
    - This will create a folder "Release" inside of the "build" folder.
    - Inside the "Release" folder, the executable "TPMS-Executable.exe" can be found.
- Linux:
    - The executables will be made directly in the build folder, there is not a "Release" folder.
- **HPC**
    - In HPC mode, please change CMakeLists.txt file. For regular computers, using the CUDA compilation flag to "all-major" doesn't seem to be a problem. In HPC, please change that to the compute capability of the available GPU. If the GPU is computer capability 7.0, then change from "all-major" to "70". 

The CMAKE file is the same for Windows and Linux. It has been tested on Windows 11, Ubuntu (22.04 or more recent), and Rocky Linux (HPC). So far it seems to work everywhere, save a few issues when running on PowerShell.

# Publications

There are currently no publications associated with this code.

# Authors

- Main developer: Andre Adam (The University of Kansas)
    - [ResearchGate](https://www.researchgate.net/profile/Andre-Adam-2)
    - [GoogleScholar](https://scholar.google.com/citations?hl=en&user=aP_rDkMAAAAJ)
    - [ORCID](https://orcid.org/0000-0002-4502-3033)
    - [GitHub](https://github.com/adama-wzr)
    - [Website](https://adama-wzr.github.io/)
- Ideator: Silven Stallard (The University of Kansas)
    - [ResearchGate](https://www.researchgate.net/profile/Silven_Stallard)
- Advisor: Dr. Guang Yang (Oak-Ridge National Laboratory)
    - [Website](https://www.ornl.gov/staff-profile/guang-yang)
    - [GoogleScholar](https://scholar.google.com/citations?user=Ph_5mDMAAAAJ&hl=en)
- Advisor: Dr. Theodore L. Bergman (The University of Kansas)
    - [GoogleScholar](https://scholar.google.com/citations?user=uHtkVqwAAAAJ&hl=en&oi=ao)
    - [ScholarGPS](https://scholargps.com/scholars/99871203759685/theodore-l-bergman)
- Advisor and Corresponding Author: Dr. Xianglin Li (Washingtion University in St. Louis)
    - [Website](https://xianglinli.wixsite.com/mysite)
    - [GoogleScholar](https://scholar.google.com/citations?user=8y0Vd8cAAAAJ&hl=en)

# Acknowledgements

This work wouldn't be possible without the computational time awarded as part of the following grants:

This work used Expanse(GPU) at SDSC through allocations MAT210014 and MAT230071 from the Advanced Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS) program, which is supported by National Science Foundation grants #2138259, #2138286, #2138307, #2137603, and #2138296.

# Upcoming Changes

- Additional capabilities are planned for the saveTPMS test function:
    - Print sub-channel information
    - Print fully-connected sub-channel information only
    - Congruency tests

# References

