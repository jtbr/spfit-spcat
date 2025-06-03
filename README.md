# SPFIT/SPCAT program for rotational spectroscopy

This repository is used for a public, version controlled set of source code for Herb M. Pickett (HMP)'s SPFIT and SPCAT programs for fitting and simulating rotational spectra. Originally called the CalPGM program suite, it also includes DPFIT and DPCAT, as well as CALMRG, although these seem no longer to be used. The original code dates to 1989.

The base version of this code was obtained from the [Cologne Database for Molecular Spectroscopy](https://cdms.astro.uni-koeln.de/classic/predictions/pickett/quelle/), with the latest modifications by Dr. Holger S. P. Müller (the `Pickett_neu.zip` file package).

This repository also includes modifications by Kelvin Lee, to allow the partition function routine to calculate the rotational partition function up to 1000 K programmatically without needing to repetitively call SPCAT.

## Build and Installation Instructions

The Makefile shows how the various files are to be linked. The programs have been tested with Microsoft Visual C++ compiler and the gnu gcc compiler (which is freely available for unix [and windows](https://cygwin.com/) platforms).

If gcc is properly configured, building should be a simple matter of running

```sh
make
```

from within the code directory. If you wish to install executables for use throughout the system, run
`make install`. To remove build products, or if you wish to force a re-build, run `make clean`.

The programs should work without modification with any ANSI compliant 'c' compiler on any size computer. All arrays are allocated dynamically, and addressing or memory limits may place a practical limit on the size of matrices that can be used. On modern computers, these should not be binding.

## Documentation and Getting Started

There is a very useful set of notes and usage information available [here](http://info.ifpan.edu.pl/~kisiel/asym/pickett/crib.htm).

You will also need to read the included [original documentation for the SPCAT and SPFIT](spinv.md).

[Original documentation for DPFIT and DPCAT](dpi.md) are also included.

## Files Contents

The identities of the included files are:

- `calfit.c`, `calcat.c`, and `calmrg.c` are the main programs.
- `subfit.c` is supplementary to calfit.
- `ulib.c`, `blas.c`, and `cnjj.c` are generic libraries.
- `calpgm.h`, `cnjj.h`, and `blas.h` are required header files.
- `slibgcc.c` contains system dependent functions.
- `spinv.c` contains functions for spins and multiple vibrations. The executables using this library and
calfit or calc at are called `spfit` and `spcat` resp e ctively.
- `dpi.c` contains functions for doublet pi with a nuclear spin The executables using this library and
calfit or calc at are called `dpfit` and `dpcat` respectively.
- `*.nam` are parameter name files for function getlbl in subfit. They are only used to label the output from calfit. The first default directory is the current directory. The second default directory is given in an environment variable named SPECNAME . Under windows, put a line like `SPECNAME=c:\spectra\` in the environment you're running in. Under Unix put a line like `SPECNAME=/home/user/spectra` in one of your initialization files (e.g. `.profile`)
- `blas.c` contains needed LINPACK double precision Basic Linear Algebra Subroutines (these may be
available on some systems in a machine coded and/or vector processor form).
- `Makefile` is the makefile for the gcc compilation.
- [`spinv.md`](spinv.md) is the specific documentation for the SPFIT and SPCAT and [`dpi.md`](dpi.md) is the specific
documentation for DPFIT and DPCAT.

## Adjustment of Marquardt-Levenberg Parameter Using a Trust-Region Approach

The *trust region* approach to Marquardt-Levenberg parameter adjustment is described in John. E. Dennis
and Robert B. Schnabel, Numerical Methods for Unconstrained Optimization and Non-linear Equations,
Prentice- Hall, 1983. The basic idea is that there is a region over which a linear least squares fit can be
*trusted*. First the value of each parameter is scaled so that the squares of the derivatives, summed over
all the lines, is unity. Then a simple least squares fit is attempted with a value of zero for the
Marquardt-Levenberg parameter, $\lambda$. If the length of the normalized parameter change vector is less than
the trust region size, then this fit is used. Otherwise, a new $\lambda$ is found in which length of the normalized
parameter change vector is equal to the trust re gion size. In the first iteration, the trust region size is set
to the length of the parameter change vector when the input value of $\lambda$ is used for the fit. For the following
iterations, the trust region is doubled if it appears that the fit is over-damped. If the fit is diverging, the
trust region is changed by a factor of 0.1 to 0.5 and the parameters from the last good fit are used again.
When the trust region is decreased, the corresponding value of $\lambda$ increases. It is a good idea to start a fit
with the Marquardt-Levenberg parameter, $\lambda$, set to zero. If the fit starts diverging, then the trust region
will be decreased appropriately.
