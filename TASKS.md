# SPFIT/SPCAT Modernization Task List

This document outlines the prioritized tasks for modernizing the SPFIT/SPCAT software suite, focusing on enabling programmatic use, improving performance, and enhancing readability/maintainability, while ensuring quantitative results and current functionality are retained.

## Completed Tasks

- [x] **Task 1: Code Readability and Clarity Improvements (Preliminary)**
    - Enhanced readability of core C code: comments, function documentation, targeted variable renames.
    - `spinv.c` split into `spinv_setup.c`, `spinv_spin_symmetry.c`, `spinv_linalg_sort.c`, `spinv_hamiltonian.c`, `spinv_utils.c` with `spinv_internal.h`.

- [x] **Task 2: Fix Compiler Warnings and Remove Unused Code**
    - All compiler warnings fixed; compiles cleanly with `-Wall -Wextra`. `cppcheck` clean. Unused functions/files removed.

- [x] **Task 3: Obsolete Functions and Portability Issues**
    - `gets()` replaced with `fgets()`. Pointer cast portability warnings fixed in `dblas.c`.

- [x] **Task 6: Build System Modernization**
    - CMake build system added alongside existing Makefile. Both are maintained independently.
    - **Critical build requirements for exact numerical reproduction** (discovered during Task 4):
        - To match v2008 test results, must use **`dblas.c` fallback** instead of OpenBLAS. OpenBLAS uses SIMD/FMA internally in `ddot`/`daxpy`, changing floating-point accumulation order at the ULP level. This flips the sign of near-zero diagonal elements in `dqrfac`'s Householder reflections, producing sign-flipped (but physically equivalent) Cholesky columns in `.var` output.
        - Must compile with **`-ffp-contract=off`** to prevent compiler FMA contraction even in simple `dblas.c` loops under `-march=native`.
        - CMake: `cmake .. -DUSE_SYSTEM_BLAS=OFF`. Makefile: leave `BLASLIB` undefined to use fallback `$(LBLAS)=dblas.o`.

- [x] **Task 4: Modularization and Decoupling with C++ Interface**
    - **All steps complete. 55/55 test files match v2008 reference baseline.**
    - Steps 1–5: `CalculationEngine` abstract interface, `SpinvEngine`/`DpiEngine` implementations with `SpinvContext`/`DpiContext` structs replacing globals, preliminary simplifications.
    - Step 6 — CalFit refactoring:
        - `CalFit` class: `initializeParameters`, `processLinesAndSetupBlocks`, `performIteration`, `finalizeOutputData`.
        - `CalFit_helpers.cpp`: `linein`, `lineix`, `getblk`, `getdbk`, `dnuadd`, `parer`, `qnfmt2`.
        - `CalFitIO`: `readInput` / `writeOutput`.
        - `fit_main.cpp` rewritten as thin wrapper.
        - `linein` accepts `std::vector<std::string>` (no `fmemopen`).
    - Step 6 — CalCat refactoring:
        - `CalCat` class: `initializeParameters`, `setupBlocks`, `computeCatalog`, `finalizeOutput`.
        - `CalCat_helpers.cpp`: `qnfmt`, `simt`, `ibufof`, `sblk_alloc`.
        - `CalCatIO`: `readInput` (reads `.int` and `.var` files).
        - `cat_main.cpp` rewritten as thin wrapper (99 lines).
        - `ibufof` block cache: static variables converted to `BlockCacheState` member struct; `tmpfile()` algorithm preserved.
    - Key bugs fixed during Step 6:
        1. **fitbgn packing mismatch**: upper-tri packed read as lower-tri. (41→46 passing)
        2. **oldpar restore missing**: unconditional `dcopy` before writing `.par`. (46→53 passing)
        3. **FMA/BLAS sign flips**: `dblas.c` + `-ffp-contract=off`. (53→55 passing)
    - **Design notes for future work**:
        - CalFit operates fully in-memory (Input/Output structs). CalCat intentionally uses `FILE*` streams for streaming output — for Python bindings (Task 9), this would need to become a callback or accumulate results in `CalCatOutput`.
        - Input structs use `std::vector`; class internals still use C-style `mallocq`/`free`, encapsulated by RAII destructors. Converting internals to `std::vector` members is possible but cosmetic — the ownership model is already correct.
        - `CalFitIO` and `CalCatIO` are independent static classes. Both call the same low-level C functions (`getpar`, `getvar`, `setopt` in `ulib.c`), but their wrapper logic differs (different file formats, second-line fields, additional data sections). Shared I/O abstraction was considered and rejected — the shared code already lives in `ulib.c`, and a higher-level wrapper would add parameterization complexity for no practical benefit.
    - **Key matrix layout facts** (documented as in-situ comments in `CalFit.cpp`):
        - `fit` is column-major, leading dimension `ndfit = nfit+1`; element `(r,c) = fit[r + c*ndfit]`.
        - `fitbgn` and `var` are packed **upper** triangular; column `j` has `j+1` elements.
        - The "copy fitbgn→fit" loop fills each **row** `n-1` of fit from column `n-1` of fitbgn.
        - After `lsqfit`, solution for parameter `k` is at `fit[nfit + k*ndfit]`.
    - Reference files: `calfit-orig.cpp` (original monolithic C, not compiled), `spinv-orig.c` (original spinv, not compiled).

- [x] **Task 7: Error Handling and Input Validation**
  - **Description**: Replace legacy slibgcc functions and add exception-based error handling for programmatic use.
  - **Completed**:
    - Exception hierarchy: `CalError → IoError / InputError / ValidationError / NumericError` (CalError.hpp).
    - `calalloc(n)` replaces `mallocq(n)`; throws `std::bad_alloc` on failure.
    - `file_helpers::open_input/open_output` replace `fopenq`; throw `IoError` on failure.
    - `file_helpers::parse_file_args` replaces `filget`.
    - `SigintFlag` RAII class replaces `rqexit(-1)`/`rqexit(0)`.
    - `CalFit::run()` and `CalCat::run()` changed from `bool` to `void`; all error paths throw.
    - `validateInput()` added to both `CalFit` and `CalCat`.
    - All 10 legacy mains (`calmrg`, `calbak`, `stark`, `sortn`, `reassign`, `termval`, `sortegy`, `iambak`, `iamcalc`, `moiam`) converted from `.c` to `.cpp` with top-level try/catch, `slibgcc` functions replaced.
    - `fit_main.cpp` and `cat_main.cpp` wrap `run()` in try/catch with structured error messages.
    - `readopt.h` given `extern "C"` guards for C++ inclusion.
    - 55/55 regression tests pass.
    - ~~Remove slibgcc.c, replacing with simple macros or standard replacements or more modern approaches.~~ Done in Task 7: `slibgcc.c` no longer compiled; all callers replaced with `file_helpers`, `calalloc`, `SigintFlag`. (from task 9.0)


## Open Tasks

- [ ] **Task X: Investigate suspected bugs in original C code**
    - These are suspicious patterns in the original Pickett code. None are triggered by existing test cases. Investigate if relevant test cases arise or before making changes to the surrounding code.
    - `calfit-orig.cpp`: "Supplied Variance" (`inpcor > 0`) path has `for (n = 1; i <= nfit; ++n)` — loop condition uses `i` but counter is `n`. Either dead code or a latent bug. The C++ refactoring in `CalFit.cpp` implements this differently. No test molecules use `inpcor > 0`.
    - `spinv_utils.c:82`: `ibtmp = idval[0]` — comment says this should probably read `idval[1]`. May be relevant for internal rotation cases.
    - `spinv_spin_symmetry.c:105`: `iis[i] = (short)ii` — possible off-by-one (`ii - 1`?).
    - `ulib.c:675`: `fabs(*pvar) < 1.01` gate on Cholesky diagonal — prevents reading correlation block if any diagonal exceeds 1.01, which shouldn't happen for a valid correlation matrix.
    - Any bugs quashed here should probably be fixed in the original (C) branch as well

- [ ] **Task 5: Performance Optimizations**
    - **Description**: Speed up calculations and exploit modern hardware.
    - **Implementation Plan**:
        - **Profile first**: Use `perf` to identify actual bottlenecks. The main cost is likely Hamiltonian computation (`hamx`/`specfc`), not BLAS. Don't optimize speculatively.
        - **Parallelization**: CalCat's outer block loop is the most promising target for OpenMP. Blocks are independent once parameters are set. CalFit's line loop accumulates into a shared matrix, making it harder to parallelize without restructuring.
        - **BLAS considerations**: `dblas.c` is required for exact v2008 reproduction (see Task 6). If profiling shows BLAS is a bottleneck, options are: (a) recompile OpenBLAS with `-ffp-contract=off -fno-associative-math`, (b) create a new baseline from v2008 code compiled with modern optimizations, or (c) accept numerically-equivalent-but-not-identical output. Profile before deciding.
        - **Memory access**: Review cache utilization in the Hamiltonian and intensity matrix computation inner loops.
    - **Constraint**: Quantitative results must remain unchanged unless a new baseline is established.


- [ ] **Task 9.0: legacy cleanup**:
    - Break apart calpgm.h header (probably into calpgm_types.h, blas_compat.h and headers for cnjj/slib/catutil/slib/ulib/subfit), including only what's needed and keeping them out of C++ headers if possible. More generally, consider usage of headers and if they should be changed or re-organized.
    - move legacy apps aside, probably into a `legacy` subdirectory, keeping only spfit/spcat as the main targets (they should remain buildable/usable). Consider other code reorganization.

- [ ] **Task 9.1: Python interface and example usage**
    - C++ example demonstrating programmatic use of CalFit/CalCat (construct engine, populate Input struct, call `run()`, read Output struct).
    - Python bindings via `pybind11` wrapping CalFit/CalCat.
    - **Design consideration**: CalCat currently streams output via `FILE*`. Python bindings will need either (a) in-memory accumulation in `CalCatOutput`, or (b) a Python callback for each catalog line. Option (a) is simpler; option (b) is more memory-efficient for large catalogs.
    - Shared I/O abstraction between CalFitIO/CalCatIO is not needed — the low-level parsing (`getpar`/`getvar` in `ulib.c`) is already shared; the wrapper differences reflect genuinely different file formats.

- [ ] **Task 10: Documentation Improvement**
    - Enhance code documentation. Existing: `spinv.md` (SPFIT/SPCAT algorithms), `dpi.md` (DPFIT/DPCAT).
    - Create user-facing documentation for the C++ API and build system.
    - Document the scientific principles and algorithms for contributors.

## Future Considerations

- [ ] **Task 11: Re-evaluate Quantum Number (QN) Limits**
    - `MAXQN` is currently 10, matching the Pickett catalog format specification. Changing it would break compatibility with existing catalogs and tools.
    - A better approach may be to support extended formats through new code paths (e.g., a wider output format option) rather than changing the constant.
    - Relevant code: `calpgm.h`, `spinv_utils.c:getqn`, `CalFit_helpers.cpp:getblk`, `CalCat.cpp:computeCatalog`.

- Any reason to consider use of modern linear algebra libraries to replace legacy code?

## Build & Test

```bash
# CMake (recommended)
cd build && cmake .. -DUSE_SYSTEM_BLAS=OFF && make spfit spcat

# Makefile (leave BLASLIB undefined for dblas.c fallback)
make spfit spcat

# Test suite
cd spfit_spcat_test_suite
python3 run_tests.py <output_subdir> && python3 compare_results.py <output_subdir>
```
