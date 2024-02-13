file(GLOB CORE_HEADERS GLOB *.hpp MovingBeam/*.hpp)
file(GLOB CORE_SOURCE *.cpp MovingBeam/*.cpp)

add_library(FinchCore ${CORE_SOURCE})
add_library(Finch::Core ALIAS FinchCore)

target_link_libraries(FinchCore Cabana::Grid nlohmann_json::nlohmann_json)

target_include_directories(FinchCore PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
install(TARGETS FinchCore DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(TARGETS FinchCore
  EXPORT Finch_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT Finch_Targets
  FILE Finch_Targets.cmake
  NAMESPACE Finch::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Finch)

install(FILES ${CORE_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
