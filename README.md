# Finch

Finite difference heat transfer using Cabana for additive manufacturing

## Dependencies

|Dependency | Version  | Required | Details|
|---------- | -------  |--------  |------- |
|CMake      | 3.16+    | Yes      | Build system
|[Kokkos](https://github.com/kokkos/kokkos) | 4.6.02 | Yes | Performance-portable execution and memory model
|[Cabana](https://github.com/ECP-copa/Cabana) | 0.8.0  | Yes | Performance portable particle/grid library
|[ADIOS2](https://github.com/ornladios/ADIOS2) | 2.11+ | No | Persistent, performance-portable BP5 field output
|[json](https://github.com/nlohmann/json)     | 3.10+   | Yes | Input files


## Build Finch

Finch requires a C++17 compiler, CMake 3.16 or newer, MPI, Kokkos 4.6.02,
Cabana 0.8.0 with `Cabana::Grid`, and nlohmann/json 3.10 or newer. On a
Debian/Ubuntu CPU system, the system prerequisites are typically
`build-essential cmake git libopenmpi-dev openmpi-bin libhwloc-dev`.

The root `Allwmake` script performs the complete reproducible build and
installs Finch and its pinned dependencies under the ignored `build/`
directory. From the repository root, run:

```bash
./Allwmake
```

The default is a Release CPU build using `mpicc`/`mpicxx`, OpenMP plus Serial
Kokkos backends, native CPU tuning, MPI-enabled Cabana Grid, static ADIOS2
BP5 support, and external nlohmann/json. It installs Finch to `build/install`.
Existing valid source and build directories are reused, so rerunning the
script performs an incremental build rather than cloning over them.

Common overrides are environment variables:

```bash
# Control parallel build jobs. WM_NCOMPPROCS is also accepted.
FINCH_JOBS=16 ./Allwmake
WM_NCOMPPROCS=16 ./Allwmake

# Select another MPI/compiler toolchain and installation location.
FINCH_CC=mpiicx FINCH_CXX=mpiicpx \
FINCH_INSTALL_PREFIX=/path/to/finch ./Allwmake

# Minimal build without ADIOS2; field output defaults to BOV/XDMF.
FINCH_ENABLE_ADIOS2=OFF ./Allwmake

# Representative CUDA and HIP builds. Use a separate dependency/build root
# when retaining more than one backend build.
FINCH_BACKEND=CUDA FINCH_CXX=/path/to/nvcc_wrapper \
FINCH_KOKKOS_ARCH=AMPERE80 \
FINCH_DEPS_ROOT="$PWD/build/dependencies-cuda" \
FINCH_BUILD_DIR="$PWD/build/cuda" ./Allwmake

FINCH_BACKEND=HIP FINCH_CXX=hipcc FINCH_KOKKOS_ARCH=VEGA90A \
FINCH_DEPS_ROOT="$PWD/build/dependencies-hip" \
FINCH_BUILD_DIR="$PWD/build/hip" ./Allwmake
```

Supported `FINCH_BACKEND` values are `OPENMP`, `SERIAL`, `CUDA`, `HIP`, and
`SYCL`. `FINCH_KOKKOS_CMAKE_ARGS`, `FINCH_CABANA_CMAKE_ARGS`,
`FINCH_ADIOS2_CMAKE_ARGS`, and `FINCH_JSON_CMAKE_ARGS` append site-specific
dependency options. Arguments passed directly to `Allwmake` are appended to
the Finch CMake configuration, for example:

```bash
./Allwmake -DCMAKE_CXX_FLAGS="-march=znver3"
```

Native CPU tuning gives the best build for the machine performing the compile;
set `FINCH_KOKKOS_ARCH=NONE` when producing a CPU binary that must run on
different architectures. For CPU builds, set `FINCH_CXX` to the MPI C++
wrapper for the desired compiler toolchain when it is not named `mpicxx`, for
example `FINCH_CXX=mpiicpx` or `FINCH_CXX=CC`.

ADIOS2 is built statically above so it can link against the default static
Kokkos installation without requiring position-independent Kokkos objects.
ADIOS2's non-CMake `adios2-config` generator is disabled because its
post-install probe does not import an external Kokkos target; Finch consumes
the supported ADIOS2 CMake package directly. `FINCH_ENABLE_ADIOS2=OFF` skips
that dependency and selects the Cabana BOV/XDMF fallback.

CMake caches both the compiler and dependency locations. `Allwmake` detects a
compiler change and refreshes only the affected generated CMake cache before
configuring. When configuring Finch manually, clear the build directory
(or select a new one) after changing toolchains:

```bash
cmake -E remove_directory build/cpu
```

For CUDA or HIP, keep the same dependency order but build Kokkos with exactly
one accelerator backend and the architecture matching the target GPU, then
use the same accelerator compiler for Cabana and Finch. Representative Kokkos
options are `Kokkos_ENABLE_CUDA=ON`, `Kokkos_ENABLE_CUDA_LAMBDA=ON`, and
`Kokkos_ARCH_AMPERE80=ON` for an A100, or `Kokkos_ENABLE_HIP=ON` and
`Kokkos_ARCH_VEGA90A=ON` for an MI200-series GPU. Do not leave
`Kokkos_ARCH_NATIVE` enabled for a cross-compiled accelerator build.

Distributed accelerator runs pass Cabana's device-resident halo buffers
directly to MPI. They therefore require a GPU-aware MPI deployment. After
verifying that support, set `FINCH_GPU_AWARE_MPI=1` for multi-rank accelerator
runs. Single-rank accelerator and ordinary CPU MPI runs do not require it.

### MPI and OpenMP placement

Use `"ranks_per_dim": [0, 0, 0]` unless a benchmark demonstrates a better
manual decomposition. Finch then chooses the Cartesian rank grid with the
smallest estimated internal halo surface for the problem geometry.

On the CPU system used for the Finch performance study, one OpenMP thread per
MPI rank was faster than hybrid layouts at both 8 and 32 total cores. A good
Open MPI starting point is therefore:

```bash
OMP_NUM_THREADS=1 OMP_PROC_BIND=close OMP_PLACES=cores \
mpirun --bind-to core --map-by core -np <ranks> \
  build/install/bin/finch -i <inputs.json>
```

Start with one rank per physical core, then sweep lower rank counts for the
actual problem. Strong scaling stops helping when local subdomains become too
small and halo exchange dominates. On multi-NUMA systems, also compare pure
MPI with one rank per NUMA domain and OpenMP threads confined to that domain.

For accelerators, use one MPI rank per GPU, one host thread per rank, and have
the scheduler or launcher bind each local rank to a distinct GPU and nearby
CPU/NUMA resources. Set `FINCH_GPU_AWARE_MPI=1` only after confirming that the
MPI library can directly communicate device buffers. Finch supports CPU and
accelerator builds from the same source, but it does not divide one
domain concurrently between an OpenMP backend and a GPU backend.

Finch labels its timestep, diffusion, heat-source, physical-boundary, halo,
runtime-function, and output regions for Kokkos Tools. Every run also reports
aggregate wall time and node-update throughput suitable for strong- and
weak-scaling comparisons.

## Run Finch

The examples and their inputs are described in the
[`examples/README`](examples/README.md). From the repository root:

```bash
examples/create_scan_path/run_example.sh
examples/single_line/run_example.sh full_physics 4
```

Generated fields and CSV files are placed under each example's ignored
`output/` directory. An ADIOS2-enabled build writes the complete temperature
time series to `fields.bp`, with embedded Fides and VTK image-data metadata for
ParaView. A build without ADIOS2 writes `grid_temperature_*.xmf` metadata and
matching raw `.dat` fields through the Cabana fallback.

## Citing

If you use Finch in your work, please cite the current release or version used from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.10698939).

## License

Finch is distributed under an [open source 3-clause BSD license](LICENSE).
