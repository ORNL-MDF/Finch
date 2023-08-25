file(GLOB HEADERS GLOB *.hpp movingBeam/*.hpp)

file(GLOB SOURCES *.cpp movingBeam/*.cpp)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(CabanaFD ${SOURCES})
target_link_libraries(CabanaFD Cabana::Cajita yaml-cpp::yaml-cpp)

target_include_directories(CabanaFD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
install(TARGETS CabanaFD DESTINATION ${CMAKE_INSTALL_LIBDIR})
