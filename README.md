# Finch

Finite difference heat transfer using Cabana for additive manufacturing

## Dependencies

|Dependency | Version  | Required | Details|
|---------- | -------  |--------  |------- |
|CMake      | 3.16+    | Yes      | Build system
|[Kokkos](https://github.com/kokkos/kokkos) | 4.6.02 | Yes | Performance-portable execution and memory model
|[Cabana](https://github.com/ECP-copa/Cabana) | 0.8.0  | Yes | Performance portable particle/grid library
|[json](https://github.com/nlohmann/json)     | 3.10+   | Yes | Input files


## Build Finch
Building Finch requires Cabana 0.8.0, Kokkos 4.6.02, MPI, a C++17 compiler, and json. A simple CPU build is shown below. Select the required Kokkos backend and architecture when building for CUDA, HIP, or SYCL; Finch uses the resulting `Kokkos::DefaultExecutionSpace` without backend-specific source changes.

```
# First build Kokkos.
# Change this path to your local Kokkos location
export KOKKOS_DIR=$HOME/kokkos
pushd $KOKKOS_DIR
mkdir build
pushd build
cmake \
    -D CMAKE_INSTALL_PREFIX=install \
    -D CMAKE_BUILD_TYPE="Release" \
    -D Kokkos_ENABLE_OPENMP=ON \
    .. ;
make -j install
popd
popd

# Next, build Cabana, pointing to the Kokkos install
# Change this path to your local Cabana location
export CABANA_DIR=$HOME/Cabana
pushd $CABANA_DIR
mkdir build
pushd build
# Note that Finch requires the Cabana::Grid sub-package, which requires MPI.
cmake \
    -D CMAKE_PREFIX_PATH=$KOKKOS_DIR/build/install \
    -D Cabana_ENABLE_GRID=ON \
    -D CMAKE_INSTALL_PREFIX=install \
    -D CMAKE_BUILD_TYPE="Release" \
    .. ;
make -j install
popd
popd

# The json library for input parsing will be automatically downloaded and included in the Finch build. 
# It can also be built externally and included in the same manner as Cabana if needed.

# Finally, build Finch
mkdir build
pushd build
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_PREFIX_PATH="$KOKKOS_DIR/build/install;$CABANA_DIR/build/install" \
  -D CMAKE_INSTALL_PREFIX=install \
  .. ;
make -j install
```

Distributed accelerator runs pass Cabana's device-resident halo buffers
directly to MPI. They therefore require a GPU-aware MPI deployment. After
verifying that support, set `FINCH_GPU_AWARE_MPI=1` for multi-rank accelerator
runs. Single-rank accelerator and ordinary CPU MPI runs do not require it.

Finch labels its timestep, diffusion, heat-source, physical-boundary, halo,
sampling, and output regions for Kokkos Tools. Every run also reports aggregate
wall time and node-update throughput suitable for strong- and weak-scaling
comparisons.

## Run Finch

The main Finch examples can be run with the scripts provided in `examples/`. Inputs are described in more detail in the [`examples/README`](examples/README.md).

This includes generating a scan path (the provided script defines path to the executable and inputs for the example):
```
cd examples/create_scan_path
./run_example.sh
```

and simulating the heat transfer for a simple path:
```
cd examples/single_line
./run_example.sh
```

## Citing

If you use Finch in your work, please cite the current release or version used from [Zenodo](https://zenodo.org/doi/10.5281/zenodo.10698939).

## License

Finch is distributed under an [open source 3-clause BSD license](LICENSE).
