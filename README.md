# LISFLOOD-FP 8.3: faster flood simulation on NVIDIA GPUs and an Apple Metal back end for the first-order finite-volume solver

LISFLOOD-FP 8.2 with faster, reproducible CUDA solvers and a new GPU back end
for Apple computers (Metal, FV1 solver). The numerical method of 8.2 is
unchanged.

Archived release: [Zenodo, DOI 10.5281/zenodo.13121119](https://doi.org/10.5281/zenodo.13121119).

## Contents

```
README.md                    this file
LICENSE                      GNU General Public License v3.0
source/                      the complete source, ready to compile
lisflood-8.2-to-8.3.patch    the same changes as one patch to the official 8.2
config.float.cmake           build settings: single precision, as the GPU solvers use
tests/                       LISFLOOD-FP's own test cases: quick check, full suite, 8.2-vs-8.3 comparison
paper/                       the cases, numbers and figures of the paper
```

`source/` is the official LISFLOOD-FP 8.2
([Zenodo 10.5281/zenodo.13121102](https://doi.org/10.5281/zenodo.13121102))
without the files the CMake build does not use (Windows and Intel library
bundles, the old makefiles, prebuilt binaries, manuals, test data and its git
history), plus the 8.3 changes. To apply the changes to a full 8.2 instead, from
its top folder: `patch -p1 -E --binary < lisflood-8.2-to-8.3.patch`. This folder
is also a git repository whose two commits are the 8.2 source and the 8.3
changes, so `git show` displays exactly what changed.

## What changes

| | in 8.2 | in 8.3 | solvers |
|---|---|---|---|
| 1 | the maximum-depth grid is kept in CPU memory, read and written over PCIe every time step, and allocated without the ghost cells the GPU writes: large grids stop with `cudaErrorIllegalAddress` | kept on the GPU at full size, copied to the CPU only when the `.max` file is written | FV1, FV2, DG2 |
| 2 | every NoData cell adds to one shared counter every time step, even when dry | dry NoData cells skip the addition | FV1, FV2, DG2, ACC |
| 3 | the flux kernels do not synchronise their threads at the end of each pass, so on grids wider or taller than 1920 cells repeated runs give different results | they do | FV1, FV2, DG2 |
| 4 | no GPU support on Apple computers | Metal back end | FV1 |
| | does not compile with CUDA 12 on Linux (an unused global `sync()` clashes with the system's) | renamed | build |
| | a stray one-line file `format` in the top folder hides the C++ header `<format>` on macOS | removed | build |
| | `generateDG2DEM` and `generateDG2start` use an uninitialised NoData value | initialised as in `lisflood` | DG2 input tools |
| | prints its internal version number, 8.1.0 | prints `LISFLOOD-FP 8.3` | start-up |

What we checked (details in the paper):

- On the test cases that come with LISFLOOD-FP, wherever 8.2 runs, 8.3 writes
  the same output files, byte for byte, for FV1, FV2, DG2 and ACC (except the
  last digit of the mass-balance log, which varies between any two runs).
- 8.2 stops with `cudaErrorIllegalAddress` under FV1, FV2 and DG2 on the
  Carlisle and EA test 5 cases (0.58 and 1.7 million cells) and on a 5.2
  million-cell grid; 8.3 runs them. Largest grids tested: 5.2 million cells
  (FV1, FV2, DG2) and 10.6 million (FV1).
- NVIDIA Compute Sanitizer finds invalid memory accesses and shared-memory
  hazards in 8.2 and none in 8.3.
- On a 24 h, 10.6 million-cell flood simulation on an NVIDIA A100, 8.3 needs
  9.4 min where 8.2 (with the maximum-depth grid allocated at full size, the
  minimum for it to run) needed 105 min. The saving depends on the grid: on a
  grid without NoData cells it is 72 % at 0.5 million cells and 85-89 % at 2 to
  8 million, and a few seconds on small grids.
- The Metal FV1 solver takes the same time steps as the CPU solver on the three
  cases of the paper and agrees with it to within single-precision rounding.
- The other GPU solvers of 8.2 (`acc_nugrid` and the adaptive MWDG2/HWFV1) are
  unchanged apart from the build fix, and were not tested.

Tested on: Red Hat Enterprise Linux 8.10, GCC 8.5, CUDA 12.6, NVIDIA A100 (40 GB);
macOS 27, Xcode 27, Apple M1 Max (32 GPU cores, 64 GB).

## Build on Linux with an NVIDIA GPU

Needs a C++14 compiler (tested with GCC 8.5), CUDA 11.4 or later (tested with 12.6),
CMake 3.13 or later, NetCDF-C, OpenMP and libnuma. On a cluster, load them as
modules, for example at ECMWF:

```bash
module load cuda/12.6 netcdf4/4.9.3 cmake/3.28.3
```

Then, from this folder, with `CMAKE_CUDA_ARCHITECTURES` set to your GPU (70
V100, 80 A100, 89 L40S, 90 H100):

```bash
cmake -S source -B build -DCMAKE_BUILD_TYPE=Release -D_CONFIG=$PWD/config.float.cmake -DCMAKE_CUDA_ARCHITECTURES=80
```

```bash
cmake --build build -j 16
```

If CMake does not find NetCDF, add
`-DNetCDF_INCLUDE_DIR=<netcdf>/include -DNetCDF_LIBRARY=<netcdf>/lib/libnetcdf.so`.

On a cluster:

- Compile on the login node, but run the check below on a GPU node, e.g.
  `srun --gres=gpu:1 --time=00:20:00 tests/quick_test.sh build/lisflood` (add
  your partition or account if the cluster needs them).
- Run `nvidia-smi` on the GPU node to see the GPU model, and make sure
  `CMAKE_CUDA_ARCHITECTURES` matches it; with the wrong number the build
  succeeds but the runs fail.
- If `lisflood` fails to start with a library error, a `LD_LIBRARY_PATH` set
  in `~/.bashrc` (e.g. for another compiler) may override the loaded modules:
  run `unset LD_LIBRARY_PATH`, then load the modules again.

## Build on a Mac with Apple Silicon

Needs Xcode with its command-line tools and, from Homebrew,
`brew install cmake libomp netcdf`:

```bash
cmake -S source -B build -DCMAKE_BUILD_TYPE=Release -D_CONFIG=$PWD/config.float.cmake -DOpenMP_CXX_FLAGS="-Xclang -fopenmp -I/opt/homebrew/opt/libomp/include" -DOpenMP_C_FLAGS="-Xclang -fopenmp -I/opt/homebrew/opt/libomp/include" -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_C_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib
```

```bash
cmake --build build -j 8 --target lisflood generateDG2DEM generateDG2start
```

The build writes `metal_fv1_kernels.metallib` next to `build/lisflood`; keep
the two together.

## Check

```bash
build/lisflood -version
```

prints `LISFLOOD-FP 8.3` among its start-up lines. Then:

```bash
tests/quick_test.sh build/lisflood
```

on an NVIDIA GPU, or `GPU=metal tests/quick_test.sh build/lisflood` on a Mac,
ends with `QUICK TEST PASSED`. See `tests/README.md` for the full suite and
`paper/README.md` to reproduce the paper.

## Use

As 8.2. Put `cuda` in the `.par` file, or pass `-cuda`, to run on the GPU; on a
Mac `metal` / `-metal` means the same. On a Mac, FV1 runs on the Apple GPU and
the other solvers on the CPU.

## Licence

GNU General Public License v3.0 (`LICENSE`), the licence cited in the
LISFLOOD-FP 8.0 and 8.2 papers and agreed with the LISFLOOD-FP developers.
