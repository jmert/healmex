# See parent directory for notes on version minimum
cmake_minimum_required(VERSION 3.9.0)
cmake_policy(SET CMP0017 NEW) # https://cmake.org/cmake/help/v3.11/policy/CMP0017.html
project(healmex CXX)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}
    "-g -O -fexceptions -fno-omit-frame-pointer")
set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS}
    "-Wl,--no-undefined")

set(MATLAB_FIND_DEBUG 1)
find_package(Matlab 9.4 REQUIRED
    COMPONENTS DATAARRAY_LIBRARY)
find_package(OpenMP)
find_package(PkgConfig REQUIRED)

# HEALPix doesn't declare its need to have cfitsio linked
pkg_search_module(HEALPIX_CXX REQUIRED IMPORTED_TARGET healpix_cxx)

set(SOURCES libhealmex.cpp)

matlab_add_mex(NAME libhealmex SHARED
    SRC ${SOURCES}
    MEX_API C++
    MATLAB_DEFAULT_RELEASE R2018a
)
set_target_properties(libhealmex PROPERTIES
    CXX_STANDARD 17
    NO_SYSTEM_FROM_IMPORTED TRUE)
target_link_libraries(libhealmex
    OpenMP::OpenMP_CXX
    PkgConfig::HEALPIX_CXX)

install(TARGETS libhealmex DESTINATION matlab)

