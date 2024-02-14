add_executable(finch SingleLayer.cpp)
target_link_libraries(finch FinchCore)
install(TARGETS finch DESTINATION ${CMAKE_INSTALL_BINDIR})
