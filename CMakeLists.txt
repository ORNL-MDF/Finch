cmake_minimum_required(VERSION 3.12)
project(Finch LANGUAGES CXX C)
set(PROJECT_VERSION "0.3.0-dev")

include(GNUInstallDirs)

find_package(Cabana 0.6.1 REQUIRED COMPONENTS Cabana::Grid)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.24")
  cmake_policy(SET CMP0135 NEW)
endif()
find_package(nlohmann_json 3.10.0 QUIET)
if(NOT NLOHMANN_JSON_FOUND)
  include(FetchContent)
  # Using most recent release here
  FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
  FetchContent_MakeAvailable(json)
  # Needs to be installed only in this case.
  install(TARGETS nlohmann_json
    EXPORT Finch_Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

find_package(Git)
if(GIT_FOUND AND IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  execute_process(
    COMMAND           ${GIT_EXECUTABLE} log --pretty=format:%H -n 1
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE   Finch_GIT_COMMIT_HASH
    )
else()
  set(Finch_GIT_COMMIT_HASH "Not a git repository.")
endif()
message(STATUS "Finch Version = '${PROJECT_VERSION}'")
message(STATUS "Finch Hash    = '${Finch_GIT_COMMIT_HASH}'")

add_subdirectory(src)
add_subdirectory(applications)
add_subdirectory(utilities)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FinchConfig.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/FinchConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/FinchConfig.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Finch)
