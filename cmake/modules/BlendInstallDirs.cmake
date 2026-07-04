include(GNUInstallDirs)

set(BLEND_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/blend" CACHE PATH "BLEND C header install directory")
set(BLEND_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "BLEND library install directory")
set(BLEND_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}" CACHE PATH "BLEND executable install directory")
set(BLEND_INSTALL_CMAKEDIR "${BLEND_INSTALL_LIBDIR}/cmake/blend" CACHE PATH "BLEND CMake package install directory")
set(BLEND_INSTALL_DOCDIR "doc" CACHE PATH "BLEND documentation install directory")
