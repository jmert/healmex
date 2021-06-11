# See parent directory for notes on version minimum
cmake_minimum_required(VERSION 3.9.0)
cmake_policy(SET CMP0017 NEW) # https://cmake.org/cmake/help/v3.11/policy/CMP0017.html
project(healmex CXX)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}
    "-g -O2 -fexceptions -fno-omit-frame-pointer")
set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS}
    "-Wl,--no-undefined -static-libstdc++")

# v3.12.0 fixes an error in dealing with default library search paths in pkg_search_module.
# (see https://gitlab.kitware.com/cmake/cmake/-/merge_requests/2027)
if (${CMAKE_VERSION} VERSION_LESS 3.12.0)
    set(FIND_PKGCONFIG CONFIGS FindPkgConfig.cmake NO_CMAKE_PATH PATHS ${CMAKE_SOURCE_DIR}/cmake/)
endif()
find_package(PkgConfig REQUIRED ${FIND_PKGCONFIG})

# v3.14.0 is the first version which becomes fully R2018a-aware with the
# C++ Data API handled. We've vendored the 3.20.0 version.
set(MATLAB_FIND_DEBUG 1)
if (${CMAKE_VERSION} VERSION_LESS 3.20.0)
    set(FIND_MATLAB CONFIGS FindMatlab.cmake NO_CMAKE_PATH PATHS ${CMAKE_SOURCE_DIR}/cmake/)
endif()
find_package(Matlab REQUIRED ${FIND_MATLAB})

find_package(OpenMP)

pkg_search_module(HEALPIX_CXX REQUIRED IMPORTED_TARGET healpix_cxx)
pkg_search_module(LIBSHARP REQUIRED IMPORTED_TARGET libsharp)

set(SOURCES libhealmex.cpp)

matlab_add_mex(NAME libhealmex SHARED
    SRC ${SOURCES}
    R2018a
)
set_target_properties(libhealmex PROPERTIES
    CXX_STANDARD 17
    NO_SYSTEM_FROM_IMPORTED TRUE)
target_link_libraries(libhealmex
    OpenMP::OpenMP_CXX
    PkgConfig::HEALPIX_CXX
    PkgConfig::LIBSHARP)

install(TARGETS libhealmex DESTINATION matlab)

