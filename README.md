# Finch

Finite difference heat transfer using Cabana for additive manufacturing

## Dependencies

|Dependency | Version  | Required | Details|
|---------- | -------  |--------  |------- |
|CMake      | 3.12+    | Yes      | Build system
|[Cabana](https://github.com/ECP-copa/Cabana) | 0.6.1  | Yes | Performance portable particle/grid library
|[json](https://github.com/nlohmann/json)     | 3.10+   | Yes | Input files


## Build Finch
Building Finch requires Cabana, Kokkos & MPI (Cabana dependencies), and json (for input files). A simple example for building on CPU is shown below.

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
  -D CMAKE_PREFIX_PATH="$CABANA_DIR/build/install" \
  -D CMAKE_INSTALL_PREFIX=install \
  .. ;
make -j install
```

## Run Finch

The main Finch examples can be run with the scripts provided in `examples/`.

This includes generating a scan path (`create_scan_path/`) and simulating the heat transfer for a simple path (`single_line/`).

## License

Finch is distributed under an [open source 3-clause BSD license](LICENSE).
