#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "BLEND::blend" for configuration "Release"
set_property(TARGET BLEND::blend APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(BLEND::blend PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libblend.2.0.0.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libblend.2.dylib"
  )

list(APPEND _cmake_import_check_targets BLEND::blend )
list(APPEND _cmake_import_check_files_for_BLEND::blend "${_IMPORT_PREFIX}/lib/libblend.2.0.0.dylib" )

# Import target "BLEND::blend_static" for configuration "Release"
set_property(TARGET BLEND::blend_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(BLEND::blend_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libblend.a"
  )

list(APPEND _cmake_import_check_targets BLEND::blend_static )
list(APPEND _cmake_import_check_files_for_BLEND::blend_static "${_IMPORT_PREFIX}/lib/libblend.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
