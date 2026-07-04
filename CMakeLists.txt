cmake_minimum_required(VERSION 3.15)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ConfigUser.cmake")
  include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/ConfigUser.cmake")
endif()
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/ConfigDefault.cmake")
set(BLEND_USER_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

project(
  ${BLEND_PACKAGE_NAME}
  VERSION ${BLEND_PACKAGE_VERSION}
  DESCRIPTION ${BLEND_PACKAGE_DESCRIPTION_SUMMARY}
  LANGUAGES C
)
set(CMAKE_INSTALL_PREFIX "${BLEND_USER_INSTALL_PREFIX}" CACHE PATH "Install path prefix" FORCE)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

include(CMakePackageConfigHelpers)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
include(BlendInstallDirs)

option(BLEND_BUILD_TESTS "Build BLEND test executables" ${BLEND_BUILD_TESTS})
option(BLEND_BUILD_SHARED "Build the BLEND shared library" ${BLEND_BUILD_SHARED})
option(BLEND_BUILD_STATIC "Build the BLEND static library" ${BLEND_BUILD_STATIC})
option(BLEND_BUILD_EXAMPLES "Build BLEND example executables" ${BLEND_BUILD_EXAMPLES})

set(BLEND_EXPORT_TARGETS)

if(BLEND_BUILD_SHARED)
  add_library(blend_shared SHARED src/blend.c)
  add_library(BLEND::blend ALIAS blend_shared)
  set_target_properties(
    blend_shared
    PROPERTIES
      OUTPUT_NAME blend
      EXPORT_NAME blend
      PUBLIC_HEADER src/blend.h
      VERSION ${BLEND_LIB_VERSION}
      SOVERSION ${BLEND_LIB_SOVERSION}
  )
  target_include_directories(
    blend_shared
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )
  target_link_libraries(blend_shared PUBLIC m)
  list(APPEND BLEND_EXPORT_TARGETS blend_shared)
endif()

if(BLEND_BUILD_STATIC)
  add_library(blend_static STATIC src/blend.c)
  if(NOT TARGET BLEND::blend)
    add_library(BLEND::blend ALIAS blend_static)
  endif()
  add_library(BLEND::blend_static ALIAS blend_static)
  set_target_properties(
    blend_static
    PROPERTIES
      OUTPUT_NAME blend
      EXPORT_NAME blend_static
      PUBLIC_HEADER src/blend.h
  )
  target_include_directories(
    blend_static
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )
  target_link_libraries(blend_static PUBLIC m)
  list(APPEND BLEND_EXPORT_TARGETS blend_static)
endif()

if(NOT BLEND_EXPORT_TARGETS)
  message(FATAL_ERROR "At least one of BLEND_BUILD_SHARED or BLEND_BUILD_STATIC must be TRUE.")
endif()

if(BLEND_BUILD_TESTS)
  enable_testing()

  add_executable(blend_api_test tests/blend_api_test.c)
  if(TARGET blend_static)
    target_link_libraries(blend_api_test PRIVATE blend_static)
  else()
    target_link_libraries(blend_api_test PRIVATE blend_shared)
  endif()
  add_test(NAME blend_api_test COMMAND blend_api_test)

endif()

if(BLEND_BUILD_EXAMPLES)
  foreach(example ex01 ex02 ex03 ex04 ex05 ex06 ex07 ex08 ex09 ex10 ex11 ex12)
    add_executable(${example} doc/examples/${example}/${example}.c)
    if(TARGET blend_static)
      target_link_libraries(${example} PRIVATE blend_static)
    else()
      target_link_libraries(${example} PRIVATE blend_shared)
    endif()
  endforeach()
endif()

install(
  TARGETS ${BLEND_EXPORT_TARGETS}
  EXPORT blend-targets
  ARCHIVE DESTINATION ${BLEND_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${BLEND_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${BLEND_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${BLEND_INSTALL_INCLUDEDIR}
)

install(
  EXPORT blend-targets
  FILE blend-targets.cmake
  NAMESPACE BLEND::
  DESTINATION ${BLEND_INSTALL_CMAKEDIR}
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/blend-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  cmake/blend-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/blend-config.cmake"
  INSTALL_DESTINATION ${BLEND_INSTALL_CMAKEDIR}
  PATH_VARS BLEND_INSTALL_INCLUDEDIR
)

configure_file(
  share/tools/blend_uninstall.sh.in
  "${CMAKE_CURRENT_BINARY_DIR}/blend_uninstall.sh"
  @ONLY
)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/blend-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/blend-config-version.cmake"
  DESTINATION ${BLEND_INSTALL_CMAKEDIR}
)

install(
  FILES ChangeLog LICENSE.TXT COPYING COPYING.LESSER
  DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/blend
)

install(
  DIRECTORY doc/examples
  DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/blend
  FILES_MATCHING
    PATTERN "*.c"
    PATTERN "*.md"
    PATTERN "*.gp"
    PATTERN "*.png"
)

install(
  PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/blend_uninstall.sh"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/tools
)
