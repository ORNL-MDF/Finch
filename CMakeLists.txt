cmake_minimum_required(VERSION 3.12)
project(Finch LANGUAGES CXX C)
set(PROJECT_VERSION "0.1.0-dev")

include(GNUInstallDirs)

find_package(Cabana REQUIRED COMPONENTS Cabana::Cajita)
find_package(yaml-cpp REQUIRED)

add_subdirectory(src)
add_subdirectory(solvers)
add_subdirectory(utilities)
