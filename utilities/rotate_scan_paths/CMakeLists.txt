cmake_minimum_required(VERSION 3.12)
project(rotate_scan_path)

include(GNUInstallDirs)

set(SOURCE_FILES rotate_scan_paths.cpp)

add_executable(rotate_scan_path ${SOURCE_FILES})

find_package(yaml-cpp REQUIRED)

target_link_libraries(rotate_scan_path yaml-cpp::yaml-cpp)

install(TARGETS rotate_scan_path DESTINATION ${CMAKE_INSTALL_BINDIR})
