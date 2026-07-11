cmake_minimum_required(VERSION 3.15)

get_filename_component(BLEND_SOURCE_DIR_REALPATH "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(BLEND_BINARY_DIR_REALPATH "${CMAKE_BINARY_DIR}" REALPATH)
if(BLEND_SOURCE_DIR_REALPATH STREQUAL BLEND_BINARY_DIR_REALPATH)
  message(FATAL_ERROR "In-source builds are not allowed. Create a separate build directory and run CMake from there.")
endif()

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

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src")
configure_file(
  src/blend_version.h.in
  "${CMAKE_CURRENT_BINARY_DIR}/src/blend_version.h"
  @ONLY
)

set(BLEND_EXPORT_TARGETS)
set(
  BLEND_PUBLIC_HEADERS
  src/blend.h
  src/blend_polygon.h
  "${CMAKE_CURRENT_BINARY_DIR}/src/blend_version.h"
)
set(
  BLEND_SOURCES
  src/blend_api.c
  src/blend_boundary.c
  src/blend_contribution.c
  src/blend_error.c
  src/blend_polygon.c
  src/blend_report.c
  src/blend_utils.c
  src/blend_window.c
)
set(
  BLEND_MODULE_SOURCES
  src/blend_docs.c
  src/blend_module.c
  src/blend_monotone.c
  src/blend_simplify.c
  src/blend_window1d.c
  src/blend_window2d.c
  src/blend_window3d.c
)

if(BLEND_BUILD_SHARED)
  add_library(blend_shared SHARED ${BLEND_SOURCES})
  add_library(BLEND::blend ALIAS blend_shared)
  set_target_properties(
    blend_shared
    PROPERTIES
      OUTPUT_NAME blend
      EXPORT_NAME blend
      PUBLIC_HEADER "${BLEND_PUBLIC_HEADERS}"
      VERSION ${BLEND_LIB_VERSION}
      SOVERSION ${BLEND_LIB_SOVERSION}
  )
  target_include_directories(
    blend_shared
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>"
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )
  target_link_libraries(blend_shared PUBLIC m)
  list(APPEND BLEND_EXPORT_TARGETS blend_shared)
endif()

if(BLEND_BUILD_STATIC)
  add_library(blend_static STATIC ${BLEND_SOURCES})
  if(NOT TARGET BLEND::blend)
    add_library(BLEND::blend ALIAS blend_static)
  endif()
  add_library(BLEND::blend_static ALIAS blend_static)
  set_target_properties(
    blend_static
    PROPERTIES
      OUTPUT_NAME blend
      EXPORT_NAME blend_static
      PUBLIC_HEADER "${BLEND_PUBLIC_HEADERS}"
  )
  target_include_directories(
    blend_static
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src>"
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  )
  target_link_libraries(blend_static PUBLIC m)
  list(APPEND BLEND_EXPORT_TARGETS blend_static)
endif()

if(NOT BLEND_EXPORT_TARGETS)
  message(FATAL_ERROR "At least one of BLEND_BUILD_SHARED or BLEND_BUILD_STATIC must be TRUE.")
endif()

add_executable(blend src/blend.c ${BLEND_MODULE_SOURCES})
target_include_directories(
  blend
  PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
if(TARGET blend_static)
  target_link_libraries(blend PRIVATE blend_static)
else()
  target_link_libraries(blend PRIVATE blend_shared)
endif()

if(BLEND_BUILD_TESTS)
  enable_testing()

  add_executable(blend_api_test tests/blend_api_test.c)
  if(TARGET blend_static)
    target_link_libraries(blend_api_test PRIVATE blend_static)
  else()
    target_link_libraries(blend_api_test PRIVATE blend_shared)
  endif()
  add_test(NAME blend_api_test COMMAND $<TARGET_FILE:blend_api_test>)
  add_test(NAME blend_cli_usage COMMAND blend)
  add_test(NAME blend_cli_help COMMAND blend --help)
  add_test(NAME blend_cli_show_modules COMMAND blend --show-modules)
  add_test(NAME blend_cli_show_windows COMMAND blend --show-windows)
  set_tests_properties(
    blend_cli_show_windows
    PROPERTIES
      PASS_REGULAR_EXPRESSION "boxcar[\r\n]+cosine[\r\n]+trapezoid"
  )
  add_test(NAME blend_cli_module_exists COMMAND blend docs =)
  add_test(
    NAME blend_window1d_usage_no_error
    COMMAND bash -c "tmp=$(mktemp -d); $<TARGET_FILE:blend> window1d > $tmp/out 2> $tmp/err; status=$?; grep -q '^blend window1d - Generate 1-D blending weights' $tmp/out && grep -q '^usage: blend window1d' $tmp/out && grep -q '^REQUIRED ARGUMENTS:' $tmp/out && grep -q '^OPTIONAL ARGUMENTS:' $tmp/out && test ! -s $tmp/err; result=$?; rm -rf $tmp; test $status -eq 0 && exit $result"
  )
  add_test(NAME blend_window1d_cosine COMMAND blend window1d -R0/4 -I1 --taper_ratio=0.25 --function=cosine)
  set_tests_properties(
    blend_window1d_cosine
    PROPERTIES
      PASS_REGULAR_EXPRESSION "0 0.345491502813[\r\n]+1 1.000000000000"
  )
  add_test(NAME blend_window1d_asymmetric COMMAND blend window1d -R0/4 -I1 --taper_ratio=0.25/0.0 --function=cosine)
  set_tests_properties(
    blend_window1d_asymmetric
    PROPERTIES
      PASS_REGULAR_EXPRESSION "0 0.345491502813[\r\n]+1 1.000000000000[\r\n]+2 1.000000000000[\r\n]+3 1.000000000000[\r\n]+4 1.000000000000"
  )
  add_test(NAME blend_window1d_arbitrary_region COMMAND blend window1d -R10/14 -I1 --taper_ratio=0.25/0.0 --function=cosine)
  set_tests_properties(
    blend_window1d_arbitrary_region
    PROPERTIES
      PASS_REGULAR_EXPRESSION "10 0.345491502813[\r\n]+11 1.000000000000[\r\n]+12 1.000000000000[\r\n]+13 1.000000000000[\r\n]+14 1.000000000000"
  )
  add_test(
    NAME blend_window1d_increment_count
    COMMAND bash -c "$<TARGET_FILE:blend> window1d -R10/14 -I0.1 --taper_ratio=0/0 | awk 'END {printf \"%d %.12g\\n\", NR, $1}'"
  )
  set_tests_properties(
    blend_window1d_increment_count
    PROPERTIES
      PASS_REGULAR_EXPRESSION "41 14"
  )
  add_test(
    NAME blend_window1d_adjust_region
    COMMAND bash -c "tmp=$(mktemp -d); $<TARGET_FILE:blend> window1d -R10/14.05 -I0.1 --taper_ratio=0/0 > $tmp/out 2> $tmp/err; status=$?; awk 'END {printf \"%d %.12g\\n\", NR, $1}' $tmp/out; grep 'adjusted region end' $tmp/err; rm -rf $tmp; exit $status"
  )
  set_tests_properties(
    blend_window1d_adjust_region
    PROPERTIES
      PASS_REGULAR_EXPRESSION "42 14.1[\r\n]+blend \\[warning\\]: window1d: adjusted region end from 14.05 to 14.1 to match increment 0.1"
  )
  add_test(
    NAME blend_window1d_query
    COMMAND bash -c "printf '0.5\\n2\\n5\\n' | $<TARGET_FILE:blend> window1d --region=0/4 --increment=1 --taper_ratio=0.25"
  )
  set_tests_properties(
    blend_window1d_query
    PROPERTIES
      PASS_REGULAR_EXPRESSION "0.5 0.672745751406[\r\n]+2 1.000000000000[\r\n]+5 0.000000000000"
  )
  add_test(
    NAME blend_window1d_blendfile_product
    COMMAND bash -c "tmp=$(mktemp -d); printf '1 4 boxcar 0/0\\n2 5 cosine 0.25/0.25\\n' > $tmp/blend.txt; $<TARGET_FILE:blend> window1d -R0/6 -I1 -B$tmp/blend.txt; status=$?; rm -rf $tmp; exit $status"
  )
  set_tests_properties(
    blend_window1d_blendfile_product
    PROPERTIES
      PASS_REGULAR_EXPRESSION "0 0.000000000000[\r\n]+1 1.000000000000[\r\n]+2 0.500000000000[\r\n]+3 1.000000000000[\r\n]+4 1.000000000000[\r\n]+5 0.500000000000[\r\n]+6 0.000000000000"
  )
  add_test(
    NAME blend_window1d_blendfile_average
    COMMAND bash -c "tmp=$(mktemp -d); printf '1 4 boxcar 0/0\\n2 5 cosine 0.25/0.25\\n' > $tmp/blend.txt; $<TARGET_FILE:blend> window1d -R0/6 -I1 -Ca -B$tmp/blend.txt; status=$?; rm -rf $tmp; exit $status"
  )
  set_tests_properties(
    blend_window1d_blendfile_average
    PROPERTIES
      PASS_REGULAR_EXPRESSION "0 0.000000000000[\r\n]+1 1.000000000000[\r\n]+2 0.750000000000[\r\n]+3 1.000000000000[\r\n]+4 1.000000000000[\r\n]+5 0.500000000000[\r\n]+6 0.000000000000"
  )
  add_test(
    NAME blend_window1d_blendfile_ignores_function_taper
    COMMAND bash -c "tmp=$(mktemp -d); printf '1 4 boxcar 0/0\\n' > $tmp/blend.txt; $<TARGET_FILE:blend> window1d -R0/5 -I1 -Ftrapezoid -T0.4 -B$tmp/blend.txt > $tmp/out 2> $tmp/err; status=$?; grep -q -- '-F/--function is ignored' $tmp/err && grep -q -- '-T/--taper_ratio is ignored' $tmp/err && grep -q '1 1.000000000000' $tmp/out; result=$?; rm -rf $tmp; test $status -eq 0 && exit $result"
  )
  add_test(NAME blend_window1d_bad_ratio COMMAND blend window1d -R0/4 -I1 --taper_ratio=0.5)
  set_tests_properties(blend_window1d_bad_ratio PROPERTIES WILL_FAIL TRUE)

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
  TARGETS blend
  RUNTIME DESTINATION ${BLEND_INSTALL_BINDIR}
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
  DESTINATION ${BLEND_INSTALL_DOCDIR}
)

install(
  DIRECTORY doc/examples
  DESTINATION ${BLEND_INSTALL_DOCDIR}
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
