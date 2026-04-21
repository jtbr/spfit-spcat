# SPFIT/SPCAT program for rotational spectroscopy

This repository is based upon source code for Herb M. Pickett's SPFIT and SPCAT programs for fitting and simulating rotational spectra. Originally called the CalPGM program suite, it also includes DPFIT and DPCAT, as well as auxiliary programs CALMRG, MOIAM, STARK, TERMVAL, SORTN, REASSIGN, SORTEGY, IAMBAK, and IAMCALC. The original code dates to 1989.

The base version of this code was obtained from the [Cologne Database for Molecular Spectroscopy](https://cdms.astro.uni-koeln.de/classic/predictions/pickett/quelle/), with the latest modifications by Dr. Holger S. P. Müller (the `Pickett_neu.zip` file package).

This repository also includes modifications by Kelvin Lee, to allow the partition function routine to calculate the rotational partition function up to 1000 K programmatically without needing to repetitively call SPCAT.

## Modernization

The codebase has been refactored from monolithic C into a modular C++ architecture while preserving exact numerical compatibility with the 2008 reference binaries across a 55-file test suite.

Key changes:
- **`CalculationEngine` interface** with `SpinvEngine` and `DpiEngine` implementations (strategy pattern), replacing global state with context structs.
- **`CalFit` class** encapsulates spfit's iterative Marquardt-Levenberg fitting. I/O separated into `CalFitIO`. Entry point `fit_main.cpp` is a thin wrapper.
- **`CalCat` class** encapsulates spcat's catalog generation. I/O separated into `CalCatIO`. Entry point `cat_main.cpp` is a thin wrapper.
- Both classes take C++ input structs and produce output structs, enabling programmatic use without files.

**Command-line usage of all executables (spfit, spcat, etc.) is unchanged from the original C code. Existing input files and workflows work without modification.**

See `TASKS.md` for the full modernization roadmap and status.


## Build and Installation Instructions

### Makefile

```sh
make            # builds spfit, spcat, calmrg (default targets)
make all        # builds all executables
make install    # installs to /usr/local/bin
make clean      # removes build products
```

### CMake

CMake provides IDE integration, dependency tracking, and configurable BLAS selection. Always run from a build subdirectory:

```sh
mkdir -p build && cd build
cmake .. -DUSE_SYSTEM_BLAS=OFF
make
```

**Do not run cmake from the project root** — it will overwrite the hand-maintained `Makefile`.

### Numerical Reproduction

To match the v2008 reference outputs exactly, the build must use the bundled `dblas.c` fallback BLAS (not OpenBLAS or other optimized BLAS) and disable FMA contraction:

- **CMake**: `cmake .. -DUSE_SYSTEM_BLAS=OFF` (default flags include `-ffp-contract=off`)
- **Makefile**: leave `BLASLIB` undefined (uses `dblas.o` automatically); `CFLAGS` includes `-ffp-contract=off`

This is necessary because FMA instructions and optimized BLAS libraries change floating-point accumulation order at the ULP level, which can flip the sign of near-zero Householder reflectors in the least-squares solver. The resulting `.var` output is physically equivalent but not bit-identical to the reference. See `TASKS.md` Task 6 for full details.

## Documentation and Getting Started

There is a very useful set of notes and usage information available [here](http://info.ifpan.edu.pl/~kisiel/asym/pickett/crib.htm).

You will also need to read the included [original documentation for SPCAT and SPFIT](spinv.md).

[Original documentation for DPFIT and DPCAT](dpi.md) are also included.

## File Contents

### Main programs

- `fit_main.cpp` — spfit entry point (thin wrapper around `CalFit`)
- `cat_main.cpp` — spcat entry point (thin wrapper around `CalCat`)
- `calmrg.c` — catalog merge utility

### Core classes (C++)

- `CalFit.cpp`, `CalFit_helpers.cpp`, `CalFit.hpp` — fitting logic (Marquardt-Levenberg iteration)
- `CalFitIO.cpp`, `CalFitIO.hpp` — spfit file I/O (`.par`, `.lin`, `.fit`, `.var`)
- `CalCat.cpp`, `CalCat_helpers.cpp`, `CalCat.hpp` — catalog generation logic
- `CalCatIO.cpp`, `CalCatIO.hpp` — spcat file I/O (`.int`, `.var`)
- `CalculationEngine.hpp` — abstract interface for Hamiltonian computation

### Calculation engines

- `SpinvEngine.cpp`, `SpinvEngine.hpp`, `SpinvContext.hpp` — engine for asymmetric/symmetric tops and general molecules
- `spinv_setup.c`, `spinv_spin_symmetry.c`, `spinv_linalg_sort.c`, `spinv_hamiltonian.c`, `spinv_utils.c`, `spinv_internal.h` — underlying C implementation (parameterized via `SpinvContext`)
- `DpiEngine.cpp`, `DpiEngine.hpp`, `DpiContext.hpp` — engine for diatomic/linear molecules with nuclear spin
- `dpi.c`, `dpi.h` — underlying C implementation (parameterized via `DpiContext`)

### Libraries

- `ulib.c` — parameter I/O (`getpar`, `getvar`, `putvar`), error analysis (`calerr`, `prcorr`)
- `lsqfit.c`, `lsqfit.h` — least-squares fitting (QR factorization, Marquardt-Levenberg solver)
- `dblas.c` — fallback BLAS routines (required for exact numerical reproduction)
- `cnjj.c`, `cnjj.h` — Clebsch-Gordan coefficients
- `slibgcc.c` — system-dependent functions
- `catutil.c`, `catutil.h` — catalog utility functions
- `subfit.c` — supplementary fitting routines
- `spinit.c`, `spinit.h` — spin initialization

### Documentation

- [`spinv.md`](spinv.md) — SPFIT/SPCAT documentation
- [`dpi.md`](dpi.md) — DPFIT/DPCAT documentation
- `TASKS.md` — modernization task list and status

### Auxiliary programs (legacy, unmodernized)

These are original Pickett utility programs. They compile and link (`make all`) but have not been modernized, tested, or verified against the current codebase. They link against the old `splib.a` static library rather than the new C++ classes.

**Data preparation / workflow:**
- `sortn.c` — standalone CLI tool for sorting catalog files by frequency (wraps the `sortn()` function from `sortsub.c`, which is also linked directly into spcat)
- `calmrg.c` — merges experimental lines into a catalog structure
- `calbak.c` — converts `.cat` format back to `.lin` for re-fitting
- `sortegy.c` — sorts energy levels by quantum number rules
- `termval.c` — associates spectral lines with upper/lower energy levels from a term list
- `reassign.c` — reassigns quantum numbers to lines via a rule file

**Internal rotation (IAM) sub-suite:**
- `moiam.c` — computes internal rotation structural parameters from molecular coordinates
- `iamcalc.c` — Hamiltonian diagonalization for internal rotation analysis
- `iambak.c` — back-transforms fitted IAM parameters between models

**Specialized:**
- `stark.c` — computes Stark effect coefficients for molecular energy levels

**Other:**
- `*.nam` — parameter name files for labeling output. Searched in the current directory, then in the directory specified by the `SPECNAME` environment variable.

## Adjustment of Marquardt-Levenberg Parameter Using a Trust-Region Approach

The *trust region* approach to Marquardt-Levenberg parameter adjustment is described in John. E. Dennis
and Robert B. Schnabel, Numerical Methods for Unconstrained Optimization and Non-linear Equations,
Prentice-Hall, 1983. The basic idea is that there is a region over which a linear least squares fit can be
*trusted*. First the value of each parameter is scaled so that the squares of the derivatives, summed over
all the lines, is unity. Then a simple least squares fit is attempted with a value of zero for the
Marquardt-Levenberg parameter, $\lambda$. If the length of the normalized parameter change vector is less than
the trust region size, then this fit is used. Otherwise, a new $\lambda$ is found in which length of the normalized
parameter change vector is equal to the trust region size. In the first iteration, the trust region size is set
to the length of the parameter change vector when the input value of $\lambda$ is used for the fit. For the following
iterations, the trust region is doubled if it appears that the fit is over-damped. If the fit is diverging, the
trust region is changed by a factor of 0.1 to 0.5 and the parameters from the last good fit are used again.
When the trust region is decreased, the corresponding value of $\lambda$ increases. It is a good idea to start a fit
with the Marquardt-Levenberg parameter, $\lambda$, set to zero. If the fit starts diverging, then the trust region
will be decreased appropriately.


### Coding Note

This modernization project began as an experiment to test coding assistants handling legacy code and porting to different languages.
It has served as a good test of the limits of current AI, particularly given the large code-base size and complexity of monolithic code.
Initially, leading LLMs and agents could often give good feedback, but struggled to complete tasks correctly, and frequently got into
coding/debugging loops from which they couldn't escape. I often needed to step in and do things correctly or interrupt and redirect.
Managing context has been and remains an issue (though it is much easier with now-smaller file sizes). But LLMs and agents have developed
impressively over the last year (Spring 2025 to Spring 2026), to the point where they are much more able to complete large and complex tasks, to include debugging well-obscured bugs, nearly without correction. It has been a fruitful experiment. Hopefully the result will also prove useful.
