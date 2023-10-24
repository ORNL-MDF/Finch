# Finch

Finite difference heat transfer using Cabana for additive manufacturing

## Dependencies

|Dependency | Version  | Required | Details|
|---------- | -------  |--------  |------- |
|CMake      | 3.20+    | Yes      | Build system
|[Cabana](https://github.com/ECP-copa/Cabana)     | master   | Yes      | Performance portable particle/grid library
|[yaml-cpp](https://github.com/jbeder/yaml-cpp)     | master   | Yes      | Input files

## Build
Building Finch requires Cabana, Kokkos (a Cabana dependency), and yaml-cpp (for input files). A simple example for building on CPU is shown below.

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
cmake \
    -D CMAKE_PREFIX_PATH=$KOKKOS_DIR/build/install \
    -D Cabana_REQUIRE_HDF5=OFF \
    -D CMAKE_INSTALL_PREFIX=install \
    -D CMAKE_BUILD_TYPE="Release" \
    .. ;
make -j install
popd
popd

# Next, build the yaml-cpp library
# Change this path to your local yaml location
export YAML_DIR=$HOME/yaml-cpp
pushd $YAML_DIR
mkdir build
pushd build
cmake \
    -D CMAKE_INSTALL_PREFIX=install \
    -D CMAKE_BUILD_TYPE="Release" \
    -D YAML_BUILD_SHARED_LIBS=off \
    .. ;
make -j install
popd
popd

# Finally, build Finch
mkdir build
pushd build
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_PREFIX_PATH="$CABANA_DIR/build/install;$YAML_DIR/build/install" \
  -D CMAKE_INSTALL_PREFIX=install \
  .. ;
make -j install
```

## License

This code is not yet openly released. It is planned to be distributed under an open source 3-clause BSD license.
