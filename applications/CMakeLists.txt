add_executable(finch SingleLayer.cpp)
target_link_libraries(finch Core)
install(TARGETS finch DESTINATION ${CMAKE_INSTALL_BINDIR})

# Optional forward-mode AD sensitivity application. Built only when the
# header-only cpp_oti_lib is available; Finch itself never depends on it.
find_path(CPP_OTI_LIB_INCLUDE_DIR otinum/otinum.hpp
  HINTS
    ${CPP_OTI_LIB_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../cpp_oti_lib/include
  DOC "Path to the cpp_oti_lib include directory")

if(CPP_OTI_LIB_INCLUDE_DIR)
  message(STATUS "Finch sensitivity: cpp_oti_lib found at ${CPP_OTI_LIB_INCLUDE_DIR}")
  add_executable(finch_sensitivity Sensitivity.cpp)
  target_link_libraries(finch_sensitivity Core)
  target_include_directories(finch_sensitivity PRIVATE
    ${CPP_OTI_LIB_INCLUDE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR})
  install(TARGETS finch_sensitivity DESTINATION ${CMAKE_INSTALL_BINDIR})
else()
  message(STATUS "Finch sensitivity: cpp_oti_lib not found, skipping finch_sensitivity")
endif()
