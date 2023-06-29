cmake_minimum_required(VERSION 3.12)
project(FD-AM)

find_package(Cabana REQUIRED COMPONENTS Cabana::Cajita)

file(GLOB HEADERS GLOB *.hpp movingBeam/*.hpp)

file(GLOB SOURCES *.cpp movingBeam/*.cpp)

include(GNUInstallDirs)
install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(CabanaFD ${SOURCES})
target_link_libraries(CabanaFD Cabana::Cajita)
target_include_directories(CabanaFD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
install(TARGETS CabanaFD DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(heat_equation heat_equation_cajita.cpp)
target_link_libraries(heat_equation CabanaFD)
install(TARGETS heat_equation DESTINATION ${CMAKE_INSTALL_BINDIR})
