# Minimum version requirements:
#   3.6.0 - IMPORTED_TARGET argument to pkg_search_modules in FindPkgConfig.cmake
#   3.8.0 - CXX_STANDARD 17 available
#   3.9.0 - imported target in FindOpenMP.cmake
cmake_minimum_required(VERSION 3.9.0)
cmake_policy(SET CMP0017 NEW) # https://cmake.org/cmake/help/v3.11/policy/CMP0017.html

project(healmat CXX)

# Can be removed once cmake_minimum_required bumped beyond 3.13 where
# R2018a-aware FindMatlab.cmake is provided.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Include environment's paths library searches
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{LIBRARY_PATH})
if (NOT WIN32)
    string(REPLACE ":" ";" CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH})
endif()

include(ExternalProject)
find_package(OpenMP)
find_package(PkgConfig REQUIRED)

pkg_search_module(CFITSIO REQUIRED cfitsio)

# Basic information required to download and compile the official HEALPix
# source code.
#
# Note that both VERSION and DATE must be updated in tandem with each release
# update.
set(HEALPIX_VERSION 3.40)
set(HEALPIX_DATE 2018Jun22)
set(HEALPIX_MD5HASH e40b4a439f34b6af11aa243751b37e1c)
# These should only need to be editted if the base URLs change.
set(HEALPIX_PKG healpix-${HEALPIX_VERSION})
set(HEALPIX_BASEURL https://managedway.dl.sourceforge.net/project/healpix/)
set(HEALPIX_FILEURL Healpix_${HEALPIX_VERSION}/Healpix_${HEALPIX_VERSION}_${HEALPIX_DATE}.tar.gz)

if (OPENMP_FOUND)
    set(HEALPIX_WITH_OPENMP "--enable-openmp")
else()
    set(HEALPIX_WITH_OPENMP "--disable-openmp")
endif()

# Because we only want the the src/cxx/autotools directory, we separate the
# build process into two pieces:
#   1. Download, unpack, and patch (if needed) the sources.
#   2. In a second step, perform the build, with the dependency configured
#      for a base path within the first target's expanded source tree.
ExternalProject_Add(healpix_download
    URL ${HEALPIX_BASEURL}${HEALPIX_FILEURL}
    URL_HASH MD5=${HEALPIX_MD5HASH}
    PREFIX ${HEALPIX_PKG}
    SOURCE_DIR "${HEALPIX_PKG}/src/healpix"
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
ExternalProject_Add(healpix
    PREFIX ${HEALPIX_PKG}
    SOURCE_DIR "${HEALPIX_PKG}/src/healpix/src/cxx/autotools"
    BUILD_IN_SOURCE YES
    CONFIGURE_COMMAND
        autoreconf -i && ./configure --prefix=<INSTALL_DIR>
            --disable-shared --with-pic ${HEALPIX_WITH_OPENMP}
    DOWNLOAD_COMMAND ""
    BUILD_COMMAND
        make
    INSTALL_COMMAND
        make install
)
add_dependencies(healpix healpix_download)
ExternalProject_Get_Property(healpix install_dir)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${install_dir})

# Since we need to have the libhealpix_cxx built and installed for the
# next step to be configurable itself, we have to run the MEX build as an
# external project as well.
ExternalProject_Add(healmex
    PREFIX "mex"
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/mex"
    CMAKE_ARGS
        "-DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}"
        "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}"
        "-DCMAKE_MODULE_PATH:PATH=${CMAKE_MODULE_PATH}"
)
add_dependencies(healmex healpix)

if (NOT (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR}))
    add_custom_command(TARGET healmex POST_BUILD
        COMMENT "Copying Matlab sources to ${CMAKE_BINARY_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/matlab/ ${CMAKE_BINARY_DIR}/matlab)
endif()
