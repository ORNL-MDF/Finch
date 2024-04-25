configure_file(Finch_Core_Config.hpp.cmakein Finch_Core_Config.hpp)
file(GLOB CORE_HEADERS GLOB *.hpp)
file(GLOB CORE_BEAM_HEADERS MovingBeam/*.hpp)
file(GLOB CORE_SOURCE *.cpp MovingBeam/*.cpp)

add_library(Core ${CORE_SOURCE})
add_library(Finch::Core ALIAS Core)

target_link_libraries(Core Cabana::Grid nlohmann_json::nlohmann_json)

target_include_directories(Core PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS Core DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(TARGETS Core
  EXPORT Finch_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT Finch_Targets
  FILE Finch_Targets.cmake
  NAMESPACE Finch::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Finch)

install(FILES ${CORE_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CORE_BEAM_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MovingBeam)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/Finch_Core_Config.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
