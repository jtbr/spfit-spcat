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
│       ├── [basename].toml           # TOML fit input (saved on first TOML run if not there)
│       ├── [basename].dipoles.toml   # TOML dipoles input (saved on first TOML run if not there)
│       ├── v2008_results/            # authoritative reference for legacy path
│       │   ├── [basename].fit
│       │   ├── [basename].cat
│       │   └── [basename].out
│       └── toml_reference/           # authoritative reference for TOML path
│           ├── [basename].fitted.toml   # spfit output: full-precision parameters + variance
│           ├── [basename].catalog.toml  # spcat output: all catalog transitions
│           ├── [basename].fit           # spfit output: human-readable fit summary
│           ├── [basename].cat           # extracted from catalog.toml by run_toml_tests.py
│           └── [basename].out           # spcat output: human-readable catalog listing
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
python3 compare_toml_outputs.py <output_subdir>                          # primary: diff TOML files
python3 compare_results.py --no-intermediates <output_subdir> toml_reference  # secondary: numeric check
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

Two complementary comparisons are run against `toml_reference/`:

**Primary — `compare_toml_outputs.py` (plain diff):**  
Diffs `mol.fitted.toml` and `mol.catalog.toml` directly.  These are the canonical
TOML outputs and contain the complete data at native double precision — parameters,
full variance matrix, all catalog transitions.  Any numerical change is caught here.

**Secondary — `compare_results.py --no-intermediates` (numeric):**  
Numerically compares the human-readable text outputs.  These cover content that has
no equivalent in the TOML files:

* `mol.fit` — the per-iteration fit log (Marquardt parameter, RMS per iteration,
  trust-region changes).  `mol.fitted.toml` captures only the final fitted state.
* `mol.out` — the energy level table, per-transition verbose listing, and partition
  function.  Not present in `mol.catalog.toml`.
* `mol.cat` — redundant with `mol.catalog.toml` but kept for `compare_results.py`
  compatibility.

`--no-intermediates` skips `.par` and `.var`, which TOML mode does not produce.

**What is not covered by either comparison:**  
Neither script cross-validates the `fitted.toml` variance matrix against the
independent v2008 `.var` baseline.  The ERR column in `.cat` (which is a quadratic
form in the variance) provides indirect evidence, but cannot catch indexing bugs in
`save_fit_output_toml` that happen to cancel in the quadratic form.
`validate_variance.py` fills this gap — see below.

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

## Parameter and variance cross-validation (`validate_variance.py`)

`validate_variance.py` cross-validates both the fitted parameter values/errors and
the variance matrix in each `toml_reference/mol.fitted.toml` against the
corresponding `v2008_results/mol.var`.

**Encoding in the `.var` file:**

| Data | Format | Precision |
|------|--------|-----------|
| Parameter value | `%21.13E` | ~14 sig figs |
| Parameter error | `%15.6E` | ~7 sig figs |
| Variance `V[i,j] / erpar[i]` | `%10.7f` | ~7 sig figs |

Variance values are written 10 characters wide, 8 per line, with no delimiter
between adjacent values.  The layout is column-major packed upper-triangular:
element `(i,j)` with `i ≤ j` is at flat index `j*(j+1)/2 + i`.

**Which parameters appear in the variance:**  
NEGBCD parameters (negative id in the `.var` file) represent component-group
continuation entries and are excluded from the variance matrix.  All positive-id
parameters — including those with small errors such as 1×10⁻³⁶ — do appear in
the variance.  In `mol.fitted.toml`, NEGBCD params appear as positive-id with
`fixed=True`.  Some `.var` files have more than one spin/symmetry header line
before the parameter list (e.g., ar-so2 has two); the script detects these by
their letter-initial first character.

**What the script checks:**

* Parameter values: `|value_toml − value_var| / |value_var| ≤ 5×10⁻¹⁴`
* Parameter errors: `|error_toml − error_var| / |error_var| ≤ 5×10⁻⁷`
* Variance: `|V_toml[i,j] / erpar[i] − stored[i,j]| ≤ 5×10⁻⁸`

```sh
python3 validate_variance.py
```

Exit code is the number of molecules with failures (0 = all pass).  All
TOML-path outputs (`.cat`, `.fit`, `.out`) have also been verified against
`v2008_results` and match except for the ERR column precision improvement
described above.

## Notes

*   Both scripts copy input files to, and run from, a temporary directory before running, so the originals are never overwritten.
*   The `v2008_results/` subdirectories are the authoritative reference for the legacy path; any change to the code must preserve them.
*   The `toml_reference/` subdirectories are the authoritative reference for the TOML path.
*   The `reference_outputs/` files were downloaded from the CDMS Pickett examples website and may differ slightly from `v2008_results/` (different compiler/era; they appear unrelated to changes made by laser_kelvin).
*   `clclo2` (general_interactions) is excluded from default runs because it is very slow.
