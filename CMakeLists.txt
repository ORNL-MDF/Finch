cmake_minimum_required(VERSION 3.12)
project(FD-AM)

find_package(Cabana REQUIRED COMPONENTS Cabana::Cajita)

add_executable(heat_equation heat_equation_cajita.cpp)

target_link_libraries(heat_equation Cabana::Cajita)

install(TARGETS heat_equation DESTINATION ${CMAKE_INSTALL_BINDIR})
