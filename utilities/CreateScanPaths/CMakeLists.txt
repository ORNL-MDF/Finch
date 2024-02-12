file(GLOB SCANPATH_HEADERS GLOB *.hpp)
file(GLOB SCANPATH_SOURCE GLOB CreateScanPaths.cpp)

add_library(ScanPaths INTERFACE)
add_library(Finch::ScanPaths ALIAS ScanPaths)

target_link_libraries(ScanPaths INTERFACE Core)

target_include_directories(Core PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS Core DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(TARGETS ScanPaths
  EXPORT Finch_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(create_scan_paths ${SCANPATH_SOURCE})
target_link_libraries(create_scan_paths ScanPaths)

install(TARGETS create_scan_paths DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${SCANPATH_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
