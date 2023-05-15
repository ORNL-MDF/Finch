export KOKKOS_INSTALL_DIR=`pwd`/kokkos/build/install
export CABANA_INSTALL_DIR=`pwd`/Cabana/build/install

mkdir build
cd build
cmake \
  -D CMAKE_BUILD_TYPE="Release" \
  -D CMAKE_PREFIX_PATH="$KOKKOS_INSTALL_DIR;$CABANA_INSTALL_DIR" \
  -D CMAKE_INSTALL_PREFIX=install \
  .. ;
make install
