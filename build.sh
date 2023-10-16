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

# Build the finite difference application
mkdir build
pushd build
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_PREFIX_PATH="$CABANA_DIR/build/install;$YAML_DIR/build/install" \
  -D CMAKE_INSTALL_PREFIX=install \
  .. ;
make -j install

