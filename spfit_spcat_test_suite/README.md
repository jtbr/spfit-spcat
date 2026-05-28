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
├── run_tests.py              # legacy-path test runner
├── run_toml_tests.py         # TOML-path test runner
├── compare_results.py        # comparison script
├── [category]/               # e.g., diatomic_molecules, asymmetric_tops
│   └── [molecule]/           # e.g., co_4, h2s
│       ├── [basename].par            # legacy SPFIT input
│       ├── [basename].lin            # legacy line input
│       ├── [basename].int            # legacy SPCAT input
│       ├── [basename].toml           # TOML fit input (saved on first TOML run)
│       ├── [basename].dipoles.toml   # TOML dipoles input (saved on first TOML run)
│       ├── v2008_results/            # authoritative reference for legacy path
│       │   ├── [basename].fit
│       │   ├── [basename].cat
│       │   └── [basename].out
│       └── toml_reference/           # authoritative reference for TOML path
│           ├── [basename].fit        # written natively by spfit in TOML mode
│           ├── [basename].cat        # extracted from [basename].catalog.toml
│           └── [basename].out        # written natively by spcat in TOML mode
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
python3 compare_results.py --no-intermediates <output_subdir> toml_reference
```

#### What the TOML path produces

Each program writes different outputs depending on mode:

| Program | File | Notes |
|---------|------|-------|
| spfit   | `mol.fit`        | Human-readable fit summary; written in both modes |
| spfit   | `mol.fitted.toml`| Full-precision parameter + variance output (TOML mode only) |
| spcat   | `mol.out`        | Human-readable catalog listing; written in both modes |
| spcat   | `mol.catalog.toml` | Catalog data including all transitions (TOML mode only) |

spcat does **not** write a separate `.cat` file in TOML mode — that data lives inside
`mol.catalog.toml`.  `run_toml_tests.py` extracts the `cat_lines` array from
`mol.catalog.toml` and writes it as `mol.cat` so that `compare_results.py` can
compare it against the reference.  The `.par`, `.var`, and `.bak` files produced
by the legacy spfit path are not written in TOML mode.

The three files compared against `toml_reference/` are therefore:

* `mol.fit` — produced natively by spfit
* `mol.out` — produced natively by spcat
* `mol.cat` — synthesized from `mol.catalog.toml` by `run_toml_tests.py`

`--no-intermediates` tells `compare_results.py` to compare only these three files,
skipping `.par` and `.var` which TOML mode does not produce.

#### TOML input files

`mol.toml` and `mol.dipoles.toml` are the permanent TOML-format inputs stored in
each molecule directory.  On the first run they are generated from the legacy
`.par`/`.lin` and `.int` files using the pickett Python library and saved for all
subsequent runs.  `run_tests.py` (legacy path) skips `.toml` files when copying
inputs to its temp directory, so the two paths remain independent.

The `toml_reference/` subdirectories are the authoritative regression baseline
for the TOML path.  They differ from `v2008_results/` in the ERR column (field 2)
of `.cat` for three molecules (h2s, o2_1and3, ch3oh); see
[Precision differences: TOML vs legacy](#precision-differences-toml-vs-legacy)
below.

## Verifying Results

`compare_results.py` compares each output file numerically against a reference
subdirectory, ignoring timestamps and using tolerances for floating-point values.

```sh
python3 compare_results.py [--no-intermediates] <output_subdir> <reference_subdir>
```

Pass `v2008_results` as the reference for legacy-path outputs, or `toml_reference`
for TOML-path outputs.  Exit code is the number of differing file pairs (0 = all pass).

## Precision differences: TOML vs legacy

The TOML path stores the variance matrix at full double precision (~17 significant
figures) in `mol.fitted.toml`.  The legacy `.var` format applies a different
treatment: `putvar` normalizes each column of the Cholesky factor by
`1/erpar[i]` (making values dimensionless, near 1.0), writes them as `%10.7f`
(7 decimal places), and `getvar` multiplies back by `erpar[i]`.  The net effect
is that the legacy path passes only ~7 significant figures of variance to `calerr`.

`calerr` computes the ERR column (field 2) of each `.cat` line as a quadratic
form in the variance matrix.  With higher-precision variance, the TOML path
produces slightly more accurate ERR values — differing from v2008 by 1 ULP in
the last printed digit for 17 lines across h2s, o2_1and3, and ch3oh.  Every
other field (frequency, log-intensity, QN columns) is bit-identical to v2008.

The `toml_reference/` outputs capture the TOML path's full-precision ERR values
and serve as the regression baseline: any change that alters `toml_reference`
outputs indicates a genuine regression in the TOML path.

To regenerate `toml_reference/` (e.g. after a code change that legitimately
changes numerical results):

```sh
uv run --project ../python python3 run_toml_tests.py toml_reference
```

## Notes

*   Both scripts copy input files to, and run from, a temporary directory before running, so the originals are never overwritten.
*   The `v2008_results/` subdirectories are the authoritative reference for the legacy path; any change to the code must preserve them.
*   The `toml_reference/` subdirectories are the authoritative reference for the TOML path.
*   The `reference_outputs/` files were downloaded from the CDMS Pickett examples website and may differ slightly from `v2008_results/` (different compiler/era; they appear unrelated to changes made by laser_kelvin).
*   `clclo2` (general_interactions) is excluded from default runs because it is very slow.
