# TPMS-Characterization
This repository is dedicated to holding code pertaining to the characterization of TPMS structures with finite thicknesses.

# Requiremets

The versions are often soft requirements, meaning they are the only versions we have tested on, but might work on older versions.

- gcc >= 11.0
- C++17 or newer
- CMAKE 3.8 or newer
- OpenMP (any recent version is fine)

**Soft Requirements**

- NVIDIA GPU (only tested with compute capability >= 7.6)
- CUDA >= 11.5 (tested all the way up to 12.8)

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

The CMAKE file is the same for Windows and Linux. It has been tested on Windows 11, Ubuntu (22.04 or more recent), and Rocky Linux (HPC). So far it seems to work everywhere, save a few issues when running on PowerShell.

# Acknowledgements

This work wouldn't be possible without the computational time awarded as part of the following grants:

This work used Expanse(GPU) at SDSC through allocations MAT210014 and MAT230071 from the Advanced Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS) program, which is supported by National Science Foundation grants #2138259, #2138286, #2138307, #2137603, and #2138296.

# Upcoming Changes

- More test functions.
- Better separation to ensure operation in computers without CUDA.
