# SPFIT/SPCAT Modernization Task List

This document outlines the prioritized tasks for modernizing the SPFIT/SPCAT software suite, focusing on enabling programmatic use, improving performance, and enhancing readability/maintainability, while ensuring quantitative results and current functionality are retained.

## Completed Tasks

- [x] **Task 1: Code Readability and Clarity Improvements (Preliminary)**
    - Enhanced readability of core C code: comments, function documentation, targeted variable renames.
    - `spinv.c` split into `spinv_setup.cpp`, `spinv_spin_symmetry.cpp`, `spinv_linalg_sort.c`, `spinv_hamiltonian.cpp`, `spinv_utils.cpp` with `spinv_internal.h`.

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
        - `CalFitIO` and `CalCatIO` are independent static classes. Both call the same low-level C functions (`getpar`, `getvar` in `ulib.c`), but their wrapper logic differs (different file formats, second-line fields, additional data sections). Shared I/O abstraction was considered and rejected — the shared code already lives in `ulib.c`, and a higher-level wrapper would add parameterization complexity for no practical benefit.
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

- [x] **Task 9.0: legacy cleanup**:
  - Reorganized all sources into `src/{spfit,spcat,engine,splib,common,legacy_apps}/`. Include path: `-Isrc`, style `"module/file.h"` everywhere.
  - Replaced `calpgm.h` with focused headers: `calpgm_types.h` (types/macros), `blas_compat.h` (CBLAS aliases), `ulib.h` (ulib + slib + pcard), `subfit.h`, `sortsub.h`, `catutil.h` (readqn/gupfmt only). Public C++ headers (`CalFit.hpp`, `CalCat.hpp`, `CalculationEngine.hpp`) no longer transitively pull in BLAS aliases or C prototypes.
  - `slib.h` merged into `ulib.h`; `pcard` moved from `catutil.c` to `ulib.c`.
  - 55/55 regression tests pass.

- [x] **Task 9.1: Python interface and example usage**
  - C++ examples in `examples/fit_example.cpp` and `examples/cat_example.cpp` demonstrate programmatic use of CalFit/CalCat.
  - Python bindings via **nanobind** (not pybind11): smaller, faster, stable-ABI wheels. Packaged under `python/` with scikit-build-core + pyproject.toml; install with `pip install ./python`.
  - `OutputSink` abstraction (`FileSink` + `MemorySink`) replaces `FILE*` in CalCat's public API, enabling in-memory catalog capture without POSIX-only `open_memstream`/`fmemopen`.
  - High-level Python API: `pickett.fit_files(base_path)` and `pickett.cat_files(base_path)`.
  - Low-level Python API: `FitSession` / `CatSession` objects with `.run()`.
  - Exception hierarchy mapped to Python: `CalError`, `IoError`, `InputError`, `ValidationError`, `NumericError`.
  - 18/18 smoke tests pass (`uv run --directory python pytest tests/`).
  - See `API.md` for full documentation.


## Open Tasks

- [ ] **Task X: Investigate suspected bugs in original C code**
    - These are suspicious patterns in the original Pickett code. None are triggered by existing test cases. Investigate if relevant test cases arise or before making changes to the surrounding code.
    - `calfit-orig.cpp`: "Supplied Variance" (`inpcor > 0`) path has `for (n = 1; i <= nfit; ++n)` — loop condition uses `i` but counter is `n`. Either dead code or a latent bug. The C++ refactoring in `CalFit.cpp` implements this differently. No test molecules use `inpcor > 0`.
    - `spinv_utils.cpp:82`: `ibtmp = idval[0]` — comment says this should probably read `idval[1]`. May be relevant for internal rotation cases.
    - `spinv_spin_symmetry.cpp:105`: `iis[i] = (short)ii` — possible off-by-one (`ii - 1`?).
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
    - Any reason to consider use of modern linear algebra libraries to replace legacy code?


- [ ] **Task 10: Documentation Improvement**
    - Enhance code documentation. Existing: `spinv.md` (SPFIT/SPCAT algorithms), `dpi.md` (DPFIT/DPCAT).
    - Create user-facing documentation for the C++ API and build system.
    - Document the scientific principles and algorithms for contributors.

- [ ] **Task 12: Typed-struct input API (file-free programmatic use)**
    - **Goal**: Define typed C++ input records as the canonical public API so that `CalFit::run` and `CalCat::run` can be driven entirely from in-memory structs — no `.par`/`.lin`/`.int`/`.var` files required. The legacy file-reading path becomes a backward-compat parser layer on top of the same core.
    - **Motivation**: Task 9.1 left `CalFitIO::readInput` and `CalCatIO::readInput` opening files with `fopen()`, calling `setopt(FILE*, …)` directly, and using `getpar`/`getvar(FILE*, …)` (`src/splib/ulib.c:576,658`). The internal structs (`CalFitInput::idpar_data`) are BCD-packed bytes, not values a Python user would naturally construct. Task 12 fixes this, and also enables a future clean file format (TOML/JSON as a serialization of the same structs).
    - **Schema** (new header `src/api/InputSchema.hpp`):
        - `Parameter { int64_t id; double value, a_priori_error; bool fixed; string label; }`
        - `LineRecord { array<int,MAXQN*2> qn; int nqn; double freq, err, weight; string blend_tag; }`
        - `DipoleMoment { int64_t id; double value; bool starts_new_component; }`
        - `VibState { int index; bool oblate; int knmin,knmax,iax,iwtpl,iwtmn,ewt0; double vsym; vector<int> nuclear_spins; }` (SPINV only)
        - `SpinvOptions { int ixz,idiag,phase_flags; vector<VibState>; string nam_file; }`
        - `DpiOptions { int isdgn, nvib; }` — only two integers
        - `EngineOptions { Kind kind; SpinvOptions spinv; DpiOptions dpi; }` (discriminated union)
        - `FitInput { title, n_iterations, marquardt_param, max_obs_calc_err, param_err_scale, freq_scale, max_lines, EngineOptions, vector<Parameter>, vector<double> variance, vector<LineRecord> }`
        - `CatInput { title, CatControl, vector<DipoleMoment>, EngineOptions, vector<Parameter>, vector<double> variance }`
        Validate against file formats in spinv.md and dpi.md.
    - **Phase 1 — Schema + builders** (`src/api/builders.{hpp,cpp}`):
        - `build_fit_input(FitInput, CalculationEngine&) → CalFitInput`; `build_cat_input(CatInput, …) → CalCatInput`.
        - New pure-virtual `CalculationEngine::apply_options(const EngineOptions&)` implemented by lifting the post-parse mutation block of `setopt` / `setopt_dpi` (`spinv_setup.cpp:612-720`, `dpi.cpp:655-679`) into helpers called from both the new method and the legacy `setopt(FILE*, …)`.
    - **Phase 2 — Parsers** (`src/api/legacy_parser.{hpp,cpp}`):
        - `parse_spinv_option_lines`, `parse_dpi_option_line`, `parse_parameter_lines`, `parse_variance_lines`, `parse_dipole_lines`, `parse_line_records` — each a pure function on `vector<string>`.
        - `CalFitIO::readInput` / `CalCatIO::readInput` rewritten as: slurp file → split lines → parsers → `FitInput`/`CatInput` → `build_*input`.
    - **Phase 3 — Python bindings**: bind `FitInput`, `CatInput`, `EngineOptions`, `Parameter`, `LineRecord`, `DipoleMoment`, `VibState`, `SpinvOptions`, `DpiOptions`; add `FitSession.from_input(FitInput)` / `CatSession.from_input(CatInput)`; keep file-path constructors (now go through parser internally).
    - **Phase 4 — Tests**: 55-molecule regression remains the gate; add `python/tests/test_typed_input.py` with round-trip (parse→build→run) and typed-only (hand-constructed structs) tests.
    - **Phase 5 — API Docs**: Update `API.md` with new input structures and constructors and `README.md` with synopsis.
    - **Phase 6 — (Optional) clean file format**: Make new version of executables that accept TOML/JSON/similar format inputs and produce same-format outputs
    - **Key files**: new `src/api/InputSchema.hpp`, `src/api/builders.{hpp,cpp}`, `src/api/legacy_parser.{hpp,cpp}`; modify `src/engine/CalculationEngine.hpp`, `SpinvEngine.{hpp,cpp}`, `DpiEngine.{hpp,cpp}`, `spinv_setup.cpp`, `dpi.cpp`, `CalFitIO.{hpp,cpp}`, `CalCatIO.{hpp,cpp}`, `CalFit_helpers.cpp`, `python/src/bindings.cpp`, `python/pickett/__init__.py`, `API.md`.
    - **Constraint**: 55-file regression baseline must remain bit-identical.

## Future Considerations

- [ ] **Task 11: Re-evaluate Quantum Number (QN) Limits**
    - `MAXQN` is currently 10, matching the Pickett catalog format specification. Changing it would break compatibility with existing catalogs and tools.
    - A better approach may be to support extended formats through new code paths (e.g., a wider output format option) rather than changing the constant.
    - Relevant code: `calpgm_types.h` (MAXQN), `spinv_utils.cpp:getqn`, `CalFit_helpers.cpp:getblk`, `CalCat.cpp:computeCatalog`.
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
