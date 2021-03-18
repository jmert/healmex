# healmex

A package providing Matlab bindings to the [HEALPix][healpix] C++ package.

## Build Requirements

- GCC 7+ (or a compiler which supports C++17's structured bindings)
- Matlab R2018a+
- CMake 3.9.0+
- `pkg-config` (or `pkgconf`)
- Internet connection (to download HEALPix sources)
  - HEALPix v3.60 is used
- HEALPix's build requirements, such as:
  - GNU Autotools (`autoreconf`, `make`, etc)
  - CFITSIO

At this time, compilation and use on Linux is the only tested and supported
platform.

It is not planned to provide support for any versions of Matlab older than
R2018a since it is the first version which provides the [interleaved][]
complex data storage format. This new option finally provides a format
compatible with standard complex data types (e.g. `double complex` in C and
`complex<double>` in C++). In Matlab's prior format, complex arrays must be
duplicated (doubling memory requirements) and copied on each entrance/exit
across the Matlab/C++ boundary.

## Installation Instructions

The repository can be cloned via `git` from GitHub:
```bash
$ git clone https://github.com/jmert/healmex.git
```
Configuration is performed via CMake's `cmake` or `ccmake` commands. For
non-interactive use, the installation prefix (where the compiled MEX function
and associated Matlab scripts) are installed to can be configured with:
```bash
$ cmake -DCMAKE_INSTALL_PREFIX=~/matlab/healmex
```
where `~/matlab/healmex` can be any other path of your choice. Then to build
and install
```bash
$ make && make install
```
The build will start by downloading HEALPix, unpacking it, and compiling the
C++ library. (The C and Fortran components are not built!) After that, the MEX
function bindings to HEALPix's C++ APIs are compiled.

To use within Matlab, the installation path should be added into Matlab's
runtime path:
```matlab
>> addpath('~/matlab/healmex')
```
Adding a similar line to your `startup.m` file may be useful.

### Configuration Options

The following are available options which can be set during configuration.
Each should be set in either the interact `ccmake` prompt or via a CMake
variable definition with `cmake -DVAR=VALUE`.

* **`ASPACKAGE`**: Defaults to `ON`. If `ON`, the package is installed as a
  Matlab package. Otherwise if `OFF`, the package is installed as a flat
  set of functions.

## Usage

If installed as a Matlab package, the bindings in Matlab are presented as
methods within the `healmex` package namespace.

Help is available for each binding
```matlab
>> help healmex.pix2ang
  [theta, phi] = pix2ang(nside, ipix, varargin)

  INPUTS
    nside       The HEALPix Nside parameter.
    ipix        Pixel indices.
  ...
```
with descriptions of the calling convention for each binding. For example,
a call to `pix2ang` to convert pixel 72 in an Nside = 4 map to θ and φ
coordinates is accomplished with a call like:
```matlab
>> [th,ph] = healmex.pix2ang(4, 74)

th =

    1.4033


ph =

    0.7854

```
If installed in the flat form, the `healmex.` prefix from calls should
be removed.

[healpix]: https://healpix.sourceforge.io/index.php
[interleaved]: https://www.mathworks.com/help/matlab/matlab_external/matlab-support-for-interleaved-complex.html

