# Minimum version currently set by: FindMatlab.cmake module
cmake_minimum_required(VERSION 3.3.0)
project(healmat CXX)

include(ExternalProject)
find_package(Matlab 9.4 REQUIRED) # 9.4 == R2018a
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
