file(GLOB ROTATE_HEADERS rotate_scan_paths.hpp)
file(GLOB ROTATE_SOURCE rotate_scan_paths.cpp)

add_library(RotateScanPaths INTERFACE)
target_link_libraries(RotateScanPaths INTERFACE Core)

target_include_directories(Core PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS Core DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(TARGETS RotateScanPaths
  EXPORT Finch_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(rotate_scan_paths ${ROTATE_SOURCE})
target_link_libraries(rotate_scan_paths RotateScanPaths)

install(TARGETS rotate_scan_paths DESTINATION ${CMAKE_INSTALL_BINDIR})
