add_library(CreateScanPaths INTERFACE)
target_link_libraries(CreateScanPaths INTERFACE Core)

target_include_directories(Core PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS Core DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(TARGETS CreateScanPaths
  EXPORT Finch_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(create_scan_paths CreateScanPaths.cpp)
target_link_libraries(create_scan_paths CreateScanPaths)

install(TARGETS create_scan_paths DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES CreateScanPaths.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
