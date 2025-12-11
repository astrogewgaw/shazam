cmake_minimum_required(VERSION 3.15...3.27)

project(
  shazam
  VERSION 1.0
  LANGUAGES CXX
  DESCRIPTION "I/O for SPOTLIGHT ring buffers, with the power of SHAZAM!"
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
include(CPM)
include(LTO)
include(Misc)
include(Colors)
include(ConfigSafeGuards)

find_lto(CXX)

CPMAddPackage("gh:TheLartians/PackageProject.cmake@1.8.0")

add_library(
  ${PROJECT_NAME}
  ${CMAKE_CURRENT_SOURCE_DIR}/include/shazam/hdr.h ${CMAKE_CURRENT_SOURCE_DIR}/include/shazam/tel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/shazam/frb.h ${CMAKE_CURRENT_SOURCE_DIR}/src/hdr.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tel.cc ${CMAKE_CURRENT_SOURCE_DIR}/src/frb.cc
)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
target_link_libraries(${PROJECT_NAME})
target_include_directories(
  ${PROJECT_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include/${PROJECT_NAME}-${PROJECT_VERSION}>
)

string(TOLOWER ${PROJECT_NAME}/version.h VERSION_HEADER_LOCATION)
packageProject(
  NAME ${PROJECT_NAME}
  NAMESPACE ${PROJECT_NAME}
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
  BINARY_DIR ${PROJECT_BINARY_DIR}
  INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include
  VERSION_HEADER "${VERSION_HEADER_LOCATION}"
  INCLUDE_DESTINATION include/${PROJECT_NAME}-${PROJECT_VERSION}
)

set(PYTHON_EXT_NAME core)
find_package(
  Python 3.10 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule
)
execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE
  OUTPUT_VARIABLE NB_DIR
)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)
nanobind_add_module(
  ${PYTHON_EXT_NAME}
  STABLE_ABI
  ${CMAKE_CURRENT_SOURCE_DIR}/include/shazam/hdr.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/shazam/tel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/include/shazam/frb.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/hdr.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/tel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/frb.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/python/shazam/bind.cc
)
target_link_libraries(${PYTHON_EXT_NAME} PUBLIC)
install(TARGETS ${PYTHON_EXT_NAME} LIBRARY DESTINATION ${PROJECT_NAME})
