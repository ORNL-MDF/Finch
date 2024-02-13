file(GLOB SCANPATH_HEADERS GLOB *.hpp)
file(GLOB SCANPATH_SOURCE GLOB Finch_CreateScanPaths.cpp)

add_library(FinchScanPaths INTERFACE)
add_library(Finch::ScanPaths ALIAS FinchScanPaths)

target_link_libraries(FinchScanPaths INTERFACE FinchCore)

target_include_directories(FinchScanPaths INTERFACE
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS FinchScanPaths DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(TARGETS FinchScanPaths
  EXPORT Finch_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_executable(create_scan_paths ${SCANPATH_SOURCE})
target_link_libraries(create_scan_paths FinchScanPaths)

install(TARGETS create_scan_paths DESTINATION ${CMAKE_INSTALL_BINDIR})

install(FILES ${SCANPATH_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
