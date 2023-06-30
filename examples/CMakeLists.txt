add_executable(cbnFD-AM main.cpp)
target_link_libraries(cbnFD-AM CabanaFD)
install(TARGETS cbnFD-AM DESTINATION ${CMAKE_INSTALL_BINDIR})
