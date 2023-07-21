cmake_minimum_required(VERSION 3.12)
project(FD-AM)

include(GNUInstallDirs)

find_package(Cabana REQUIRED COMPONENTS Cabana::Cajita)

add_subdirectory(src)
add_subdirectory(examples)
