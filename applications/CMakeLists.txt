add_executable(finch SingleLayer.cpp)
target_link_libraries(finch Core)
install(TARGETS finch DESTINATION ${CMAKE_INSTALL_BINDIR})

# Optional forward-mode AD sensitivity application. Sparrow is header-only, so
# its include directory is the complete dependency integration; no Sparrow
# library target is built or linked, and Finch itself never depends on it.
find_path(SPARROW_INCLUDE_DIR otinum/otinum.hpp
  HINTS
    ${SPARROW_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Sparrow/include
  DOC "Path to the Sparrow include directory")

if(SPARROW_INCLUDE_DIR)
  message(STATUS "Finch sensitivity: Sparrow found at ${SPARROW_INCLUDE_DIR}")
  add_executable(finch_sensitivity SingleLayer_Sensitivity.cpp)
  target_link_libraries(finch_sensitivity Core)
  target_compile_definitions(finch_sensitivity PRIVATE OTI_ENABLE_KOKKOS)
  target_include_directories(finch_sensitivity PRIVATE
    ${SPARROW_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/integrations/sparrow)
  install(TARGETS finch_sensitivity DESTINATION ${CMAKE_INSTALL_BINDIR})
else()
  message(STATUS "Finch sensitivity: Sparrow not found, skipping finch_sensitivity")
endif()
