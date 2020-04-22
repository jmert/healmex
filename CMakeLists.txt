# Minimum version requirements:
#   3.6.0 - IMPORTED_TARGET argument to pkg_search_modules in FindPkgConfig.cmake
#   3.8.0 - CXX_STANDARD 17 available
#   3.9.0 - imported target in FindOpenMP.cmake
#   3.12.0 - Do not ignore default library search paths in pkg_search_module
#            (see https://gitlab.kitware.com/cmake/cmake/-/merge_requests/2027)
cmake_minimum_required(VERSION 3.12.0)
cmake_policy(SET CMP0017 NEW) # https://cmake.org/cmake/help/v3.11/policy/CMP0017.html

project(healmex CXX)
option(ASCLASS "install healmex as a Matlab class" ON)

# cmake/FindMatlab.cmake can be deleted and this line removed once
# cmake_minimum_required bumped beyond 3.14.0 where FindMatlab.cmake
# becomes fully R2018a-aware and C++ Data API compatible.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Include environment's paths library searches
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} $ENV{LIBRARY_PATH})
if (NOT WIN32)
    string(REPLACE ":" ";" CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH}")
endif()

include(ExternalProject)
find_package(OpenMP)
find_package(PkgConfig REQUIRED)

# The HEALPix package does its own search, but we can do a pkg-config search
# here during config to error out earlier.
pkg_search_module(CFITSIO REQUIRED cfitsio)

# Basic information required to download and compile the official HEALPix
# source code.
#
# Note that both VERSION and DATE must be updated in tandem with each release
# update.
set(HEALPIX_VERSION 3.60)
set(HEALPIX_DATE 2019Dec18)
set(HEALPIX_MD5HASH 9b51b2fc919f4e70076d296826eebee0)
# These should only need to be editted if the base URLs change.
set(HEALPIX_PKG healpix-${HEALPIX_VERSION})
set(HEALPIX_BASEURL https://managedway.dl.sourceforge.net/project/healpix/)
set(HEALPIX_FILEURL Healpix_${HEALPIX_VERSION}/Healpix_${HEALPIX_VERSION}_${HEALPIX_DATE}.tar.gz)

if (OPENMP_FOUND)
    set(HEALPIX_WITH_OPENMP "--enable-openmp")
else()
    set(HEALPIX_WITH_OPENMP "--disable-openmp")
endif()

# Because we only want to build the src/cxx/autotools directory, we separate
# the build process into two pieces:
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
ExternalProject_Add(libsharp
    PREFIX ${HEALPIX_PKG}
    SOURCE_DIR "${HEALPIX_PKG}/src/healpix/src/common_libraries/libsharp"
    BUILD_IN_SOURCE YES
    CONFIGURE_COMMAND
        ./configure --prefix=<INSTALL_DIR> --disable-shared --with-pic ${HEALPIX_WITH_OPENMP}
    DOWNLOAD_COMMAND ""
    BUILD_COMMAND
        $(MAKE)
    INSTALL_COMMAND
        $(MAKE) install
)
ExternalProject_Add(healpix
    PREFIX ${HEALPIX_PKG}
    SOURCE_DIR "${HEALPIX_PKG}/src/healpix/src/cxx"
    BUILD_IN_SOURCE YES
    CONFIGURE_COMMAND
        ${CMAKE_COMMAND} -E env "PKG_CONFIG_PATH=$ENV{PKG_CONFIG_PATH}:<INSTALL_DIR>/lib/pkgconfig"
        ./configure --prefix=<INSTALL_DIR> --disable-shared --with-pic ${HEALPIX_WITH_OPENMP}
    DOWNLOAD_COMMAND ""
    BUILD_COMMAND
        $(MAKE)
    INSTALL_COMMAND
        $(MAKE) install
)
add_dependencies(libsharp healpix_download)
add_dependencies(healpix libsharp healpix_download)
ExternalProject_Get_Property(libsharp install_dir)
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
add_dependencies(healmex libsharp healpix)

if (ASCLASS)
    set(PKGBASE "")
    set(CLASSPREFIX "healmex.")
else()
    set(PKGBASE "@healmex/")
    set(CLASSPREFIX "")
endif()

file(GLOB MFILES  "${CMAKE_SOURCE_DIR}/matlab/@healmex/*.m")
file(GLOB INFILES "${CMAKE_SOURCE_DIR}/matlab/@healmex/*.m.in")

# Copy Matlab scripts into the build directory.
foreach(INFILE ${INFILES})
    get_filename_component(FILENAME ${INFILE} NAME_WE)
    set(FILENAME "${CMAKE_BINARY_DIR}/matlab/@healmex/${FILENAME}.m")
    message(STATUS "Configuring: ${FILENAME}")
    configure_file(${INFILE} ${FILENAME} @ONLY)
endforeach()
if (NOT (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR}))
    foreach(MFILE ${MFILES})
        get_filename_component(FILENAME ${MFILE} NAME)
        set(FILENAME "${CMAKE_BINARY_DIR}/matlab/@healmex/${FILENAME}")
        message(STATUS "Staging: ${FILENAME}")
        execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${MFILE} ${FILENAME})
    endforeach()
endif()

# Installation instructions for MEX + scripts
install(DIRECTORY ${CMAKE_BINARY_DIR}/matlab/
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING
    PATTERN "@healmex" EXCLUDE
    PATTERN "*.mex*"
    )
install(DIRECTORY ${CMAKE_BINARY_DIR}/matlab/${PKGBASE}
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.m"
    )
