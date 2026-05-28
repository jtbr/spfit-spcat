# SPFIT/SPCAT Test Suite

This directory contains a collection of test cases for the SPFIT and SPCAT programs,
based on examples from the [CDMS Pickett software examples page](https://cdms.astro.uni-koeln.de/classic/predictions/pickett/beispiele/).

## Purpose

The purpose of this suite is to:
1. Provide a standardized set of inputs for testing SPFIT and SPCAT installations.
2. Allow verification of results by comparing generated outputs against reference outputs.

## Directory Structure

The test suite is organized as follows:

```
spfit_spcat_test_suite/
├── README.md
├── run_tests.py          # legacy-path test runner
├── run_toml_tests.py     # TOML-path test runner
├── compare_results.py    # comparison script
├── [category]/           # e.g., diatomic_molecules, asymmetric_tops
│   └── [molecule]/       # e.g., co_4, h2s
│       ├── [basename].par        # legacy SPFIT input
│       ├── [basename].lin        # legacy line input
│       ├── [basename].int        # legacy SPCAT input
│       └── v2008_results/        # authoritative reference outputs
│           ├── [basename].fit
│           ├── [basename].cat
│           └── [basename].out
```

## Running Tests

**Prerequisites:**
* SPFIT and SPCAT programs pre-built in the parent directory (`../spfit`, `../spcat`)

Navigate to the `spfit_spcat_test_suite` subdirectory before running any scripts.

### Legacy path (`.par`/`.lin`/`.int` files)

```sh
cd spfit_spcat_test_suite
python3 run_tests.py <output_subdir>
python3 compare_results.py <output_subdir> v2008_results
```

Runs `spfit` and `spcat` on copies of each molecule's legacy input files and stores the outputs in `<output_subdir>/` within each molecule directory.  `clclo2` is skipped by default (very slow); pass `--all` before `<output_subdir>` to include it.

### TOML path (`mol.toml` / `mol.dipoles.toml` files)

Requires the `pickett` Python package.  Build and install it with uv from the
`python/` directory:

```sh
cd ../python && uv pip install -e .
```

Then from `spfit_spcat_test_suite/`:

```sh
uv run --project ../python python3 run_toml_tests.py <output_subdir>
python3 compare_results.py --no-intermediates <output_subdir> v2008_results
```

`run_toml_tests.py` converts each molecule's legacy inputs to TOML format using
the pickett API, runs `spfit`/`spcat` in TOML auto-detect mode, and extracts
catalog lines from the resulting `mol.catalog.toml` for comparison.

**Note:** In TOML mode `spcat` does not write a separate `.cat` file — the
`run_toml_tests.py` script extracts `cat_lines` from `mol.catalog.toml` and
writes them as `mol.cat` so `compare_results.py` can compare them against the
reference.

`--no-intermediates` tells `compare_results.py` to compare only `.fit`, `.cat`,
and `.out` (TOML mode does not produce `.par`, `.var`, or `.bak`).

## Verifying Results

`compare_results.py` compares each output file numerically against the `v2008_results` reference, ignoring timestamps and using tolerances for floating-point values.

```sh
python compare_results.py [--no-intermediates] <output_subdir> [reference_subdir]
```

Exit code is the number of differing file pairs (0 = all pass).

## Notes

*   Both scripts copy input files to, and run from, a temporary directory before running, so the originals are never overwritten.
*   The `v2008_results/` subdirectories are the authoritative reference; any change to the code must preserve them.
*   The `reference_outputs/` files were downloaded from the CDMS Pickett examples website and may differ slightly from `v2008_results/` (different compiler/era; they appear unrelated to changes made by laser_kelvin).
*   `clclo2` (general_interactions) is excluded from default runs because it is very slow.
