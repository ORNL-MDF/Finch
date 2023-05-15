cmake_minimum_required(VERSION 3.12)
project(FD-AM)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ftree-vectorize")

# Find MPI
find_package(MPI REQUIRED)
find_package(Kokkos 3.2 REQUIRED)
find_package(Cabana REQUIRED COMPONENTS Cabana::Cajita)
include_directories(${MPI_INCLUDE_PATH})

# Add the executable
add_executable(heat_equation heat_equation_cajita.cpp)

# Link MPI libraries
target_link_libraries(heat_equation ${MPI_CXX_LIBRARIES} Kokkos::kokkos Cabana::Cajita)

install(TARGETS heat_equation DESTINATION ${CMAKE_INSTALL_BINDIR})
