# Compilation
LISFLOOD-FP can be compiled on Windows or Linux using [CMake](https://cmake.org) version 3.13 or above.

## Compiling on Windows with the MSVC compiler
Launch Visual Studio 2019 and open `lisflood-fp` as a local folder.
Choose the `msvc-x64-Debug` or `msvc-x64-Release` configuration, then `Rebuild All`.

## Compiling on Windows with the Intel compiler
Close Visual Studio, then run `launch_vs2019_intel64.bat`, which assumes that Visual Studio 2019 Community Edition is installed (the `.bat` script sets the necessary environment variables for Visual Studio to locate the Intel compiler).
Choose the `intel-x64-Debug` or `intel-x64-Release` configuration, then `Rebuild All`.
Once launched from the `.bat` file, Visual Studio can compile with either MSVC or Intel compiler by choosing the appropriate configuration.

## Running or debugging on Windows
To run `lisflood.exe` from Visual Studio, switch the `Solution Explorer` to `CMake Targets View`, the `Add Debug Configuration` to `lisflood (executable)` (see the [Visual Studio documentation](https://docs.microsoft.com/en-us/cpp/build/configure-cmake-debugging-sessions) for more details).

## Compiling on Linux

Ensure a recent version of CMake is installed.
Also ensure `libnuma-dev` is installed and, for NetCDF output and dynamic rainfall support, ensure `libnetcdf-dev` is also installed.  On Ubuntu:

````bash
sudo snap install cmake --classic
sudo apt install libnuma-dev libnetcdf-dev
````

Then open a terminal at in the root `lisflood-fp` directory:

````bash
cmake -S . -B build
cmake --build build
````

The `lisflood` executable is written to the `build` directory.
If `libnuma` is installed in a non-standard location, use

````bash
cmake -S . -B build -DNUMA_ROOT=<path>
cmake --build build
````

## Customising the build
The default build configuration is given in `config.default.cmake`. To customise the build, copy and modify this file.
Then, in Windows, edit the appropriate configuration in `CMakeSettings.json`:
````json
"configurations": [
  {
    "name": "msvc-x64-Debug",
    …
    "cmakeCommandArgs": "-D_CONFIG=<filename>",
    …
  }
  …
]
````

Or, in Linux:

````bash
cmake -S . -B build -D_CONFIG=<filename>
cmake --build build
````

## NVIDIA CUDA support
CMake automatically compiles LISFLOOD-FP FV1 and DG2 CUDA solvers if the [NVIDIA CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) is installed.
To customise the CUDA compute capabilities, see `config.default.cmake`.
