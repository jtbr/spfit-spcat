# The CALPGM Program Suite


## Installation Instructions

The Makefile shows how the various files are to be linked. The programs have been tested with Microsoft Visual C++ compiler and the gnu gcc compiler (which is freely available for unix [and windows](https://cygwin.com/) platforms). The programs should work without modification with any ANSI compliant 'c' compiler on any size computer. All arrays are allocated dynamically, and addressing or memory limits will place a practical limit on the size of matrices that can be used. For 16- bit computers the address limit is equivalent to a 90 X 90 double precision matrix, while for a 32- bit computer the addressing limit is 23170 X 23170. The program has been used on a 64- bit DEC alpha computer where the addressing limit is correspondingly
larger. For both 32-bit and 64-bit computers, a more significant practical limit is usually given by the
amount of memory or the amount of disk spac e available for virtual memory.

The identities of the files are:

 - `calfit.c`, `calcat.c`, and `calmrg.c` are the main programs.
 - `subfit.c` is supplementary to calfit.
 - `ulib.c`, `blas.c`, and `cnjj.c` are generic libraries.
 - `calpgm.h`, `cnjj.h`, and `blas.h` are required header files.
 - `slibgcc.c` contains system dependent functions.
 - `spinv.c` contains functions for spins and multiple vibrations. The executables using this library and
calfit or calc at are called `spfit` and `spcat` resp e ctively.
 - `dpi.c` contains functions for doublet pi with a nuclear spin The executables using this library and
calfit or calc at are called `dpfit` and `dpcat` respectively.
 - `*.nam` are parameter name files for function getlbl in subfit. They are only used to label the
output from calfit. The first default directory is the current directory. The second default directory is given in an environment variable named SPECNAME . Under windows, put a line like `SPECNAME=c:\spectra\` in the environment you're running in. Under Unix put a line like `SPECNAME=/home/user/spectra` in one of your initialization files (e.g. `.profile`)
 - `blas.c` contains needed LINPACK double precision Basic Linear Algebra Subroutines (these may be
available on some systems in a machine coded and/or vector processor form).
 - `Makefile` is the makefile for the gcc compilation.
 - `spinv.pdf` is the specific documentation for the SPFIT and SPCAT and `dpi.pdf` is the specific
documentation for DPFIT and DPCAT.

## Adjustment of Marquardt-Levenberg Parameter Using a Trust-Region Approach

The *trust region* approach to Marquardt-Levenb erg parameter adjustment is described in John. E. Dennis
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
