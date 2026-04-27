# SPFIT/SPCAT programs for rotational spectroscopy

This repository is based upon source code for Herb M. Pickett's **SPFIT** and **SPCAT** programs for fitting and simulating rotational spectra ([H. M. Pickett, "The Fitting and Prediction of Vibration-Rotation Spectra with Spin Interactions," J. Mol. Spectrosc. 148, 371-377 (1991)](https://www.sciencedirect.com/science/article/abs/pii/002228529190393O)). Originally called the **CalPGM program suite**, it also includes **DPFIT** and **DPCAT**, as well as auxiliary programs CALMRG, MOIAM, STARK, TERMVAL, SORTN, REASSIGN, SORTEGY, IAMBAK, and IAMCALC. The original code dates to 1989.

The base version of this code was obtained from the [Cologne Database for Molecular Spectroscopy](https://cdms.astro.uni-koeln.de/classic/predictions/pickett/quelle/), with the latest modifications by Dr. Holger S. P. Müller (the `Pickett_neu.zip` file package).

This repository also includes [modifications by Kelvin Lee](https://github.com/laserkelvin/Pickett), to allow the partition function routine to calculate the rotational partition function up to 1000 K programmatically without needing to repetitively call SPCAT.

## Modernization

The codebase has been refactored from monolithic C into a modular C++ architecture while preserving exact numerical compatibility with the 2008 reference binaries (as modified by Kelvin Lee) across a 55-file test suite.

Key changes:

- **`CalculationEngine` interface** with `SpinvEngine` and `DpiEngine` implementations (strategy pattern), replacing global state with context structs.
- **`CalFit` class** encapsulates spfit's iterative Marquardt-Levenberg fitting. I/O separated into `CalFitIO`. Entry point `fit_main.cpp` is a thin wrapper.
- **`CalCat` class** encapsulates spcat's catalog generation. I/O separated into `CalCatIO`. Entry point `cat_main.cpp` is a thin wrapper.
- Both classes take C++ input structs and produce output structs, enabling programmatic use without files.

**Command-line usage of all executables (spfit, spcat, etc.) is unchanged from the original C code. Existing input files and workflows work without modification.** New or improved interfaces may be considered for the future however.

See [TASKS.md](TASKS.md) for the full modernization roadmap and status.


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

To get started, you'll need to reference the full documentation explaining inputs and outputs for each of the programs:
 - [SPCAT and SPFIT](spinv.md)
 - [DPFIT and DPCAT](dpi.md)

These are conversions of the original PDF-format documentation into markdown format; this was not trivial so in case of doubt, refer to the originals and please report any mistakes.

There is a modern [python wrapper](https://github.com/Ltotheois/Pyckett) available for the original tools that you might find useful and should work the same with this version.


## Source Layout

Sources are organized under `src/`:

```
src/
├── spfit/       fit_main.cpp, CalFit.{cpp,hpp}, CalFit_helpers.cpp,
│                CalFitIO.{cpp,hpp}, subfit.{cpp,h}
├── spcat/       cat_main.cpp, CalCat.{cpp,hpp}, CalCat_helpers.cpp,
│                CalCatIO.{cpp,hpp}, sortsub.{c,h}
├── engine/      CalculationEngine.hpp,
│                SpinvEngine.{cpp,hpp}, SpinvContext.hpp,
│                spinv_setup.cpp, spinv_spin_symmetry.cpp,
│                spinv_linalg_sort.c, spinv_hamiltonian.cpp,
│                spinv_utils.cpp, spinv_internal.h,
│                DpiEngine.{cpp,hpp}, DpiContext.hpp, dpi.{cpp,h},
│                spinit.{cpp,h}
├── splib/       ulib.{c,h}, cnjj.{c,h}, catutil.{c,h},
│                lsqfit.{c,h}, cblas.h, blas_compat.h, dblas.c,
│                ftran.c, calpgm_types.h
├── common/      CalError.hpp, file_helpers.{cpp,hpp},
│                Logger.{cpp,hpp}, SigintFlag.{cpp,hpp}
└── legacy_apps/ calmrg.cpp, calbak.cpp, sortn.cpp, reassign.cpp,
                 sortegy.cpp, termval.cpp, stark.cpp, moiam.cpp,
                 iamcalc.cpp, iambak.cpp
```

**Key files:**

- `src/spfit/CalFit.cpp` — Marquardt-Levenberg fitting logic
- `src/spcat/CalCat.cpp` — catalog generation logic
- `src/engine/CalculationEngine.hpp` — abstract interface for Hamiltonian computation
- `src/splib/ulib.c` — parameter I/O (`getpar`, `getvar`, `putvar`), line and BCD parsing, utility functions
- `src/splib/lsqfit.c` — least-squares solver (QR factorization, Marquardt-Levenberg)
- `src/splib/cnjj.c` - Clebsch-Gordan coefficients
- `src/splib/dblas.c` — fallback BLAS (required for exact numerical reproduction of baseline)

**Runtime support files (not included in repo):**

- `*.nam` — parameter name files used by `getlbl` in `subfit.cpp` to label `.fit` output (e.g. `sping.nam`, `dpi.nam`). Searched in the current directory, then in the directory specified by the `SPECNAME` environment variable. Not included in this repository; their absence is handled gracefully (output will lack parameter name labels).

### Auxiliary programs

Legacy utility programs live in `src/legacy_apps/`. They compile and link (`make all`) but have not been modernized and the current test suite does not test them.
**Data preparation / workflow:**

- `sortn` — standalone CLI tool for sorting catalog files by frequency (wraps the `sortn()` function from `sortsub.c`, which is also linked directly into spcat)
- `calmrg` — merges experimental lines into a catalog structure
- `calbak` — converts `.cat` format back to `.lin` for re-fitting
- `sortegy` — sorts energy levels by quantum number rules
- `termval` — associates spectral lines with upper/lower energy levels from a term list
- `reassign` — reassigns quantum numbers to lines via a rule file

**Internal rotation (IAM) sub-suite:**

- `moiam` — computes internal rotation structural parameters from molecular coordinates
- `iamcalc` — Hamiltonian diagonalization for internal rotation analysis
- `iambak` — back-transforms fitted IAM parameters between models
- `ftran` — Fourier transform routines used by moiam and iamcalc
- `readopt`, `readopt.h` — option/parameter reading used by iamcalc and iambak

**Specialized:**

- `stark` — computes Stark effect coefficients for molecular energy levels

## Contributions

Contributions to this project are welcome. If you find any problems, please report them in the form of an Issue, or submit fixes as a Pull Request.

We claim no rights to the original code.


## Coding Note

This modernization project began as an experiment to test coding assistants handling legacy code and porting to different languages.
It has served as a good test of the limits of current AI, particularly given the large code-base size and complexity of monolithic code.
Initially, leading LLMs and agents could often give reasonable feedback, but struggled to complete tasks correctly, often gave poor suggestions or changed things unnecessarily and frequently got into coding/debugging loops from which they couldn't escape. I often needed to step in and do things correctly or interrupt and redirect.
Managing context has been and remains an issue (though it is much easier with now-smaller file sizes). But LLMs and agents have developed
impressively over the last year (Spring 2025 to Spring 2026), to the point where they are much more able to correctly complete large and complex tasks, to include debugging well-obscured bugs, nearly without correction. It has been a fruitful experiment. Hopefully the result will also prove useful.
