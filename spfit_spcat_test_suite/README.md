# SPFIT/SPCAT Test Suite

This directory contains regression tests for the SPFIT and SPCAT programs, based
on examples from the [CDMS Pickett software examples page](https://cdms.astro.uni-koeln.de/classic/predictions/pickett/beispiele/).

## Directory Structure

```
spfit_spcat_test_suite/
├── README.md
├── run_tests.py              # legacy-path test runner
├── run_toml_tests.py         # TOML-path test runner
├── compare_results.py        # numeric comparison (both paths)
├── compare_toml_outputs.py   # diff-based comparison (TOML path)
├── validate_variance.py      # cross-validates fitted.toml against v2008 .var
├── [category]/               # e.g., diatomic_molecules, asymmetric_tops
│   └── [molecule]/           # e.g., co_4, h2s
│       ├── [basename].par            # legacy SPFIT input
│       ├── [basename].lin            # legacy line input
│       ├── [basename].int            # legacy SPCAT input
│       ├── [basename].toml           # TOML fit input (generated on first TOML run)
│       ├── [basename].dipoles.toml   # TOML dipoles input (generated on first TOML run)
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

**Prerequisites:** SPFIT and SPCAT pre-built in the parent directory (`../spfit`,
`../spcat`, as produced by running `make`).  Run all scripts from the `spfit_spcat_test_suite/` directory.

---

## Legacy path — best practice

The legacy path runs `spfit`/`spcat` on `.par`/`.lin`/`.int` files and compares
against `v2008_results/` (the authoritative 2008 reference).

```sh
# 1. Run
python3 run_tests.py <output_subdir>

# 2. Compare against v2008 reference
python3 compare_results.py <output_subdir> v2008_results
```

`clclo2` is excluded by default (very slow); pass `--all` before `<output_subdir>`
to include it.  Exit code of `compare_results.py` is the number of differing file
pairs (0 = all pass).

---

## TOML path — best practice

The TOML path requires the `pickett` Python package (build with
`cd ../python && uv pip install -e .`).

```sh
# 1. Run spfit and spcat in TOML mode
uv run --project ../python python3 run_toml_tests.py <output_subdir>

# 2. Primary regression check: diff the TOML output files against toml_reference
python3 compare_toml_outputs.py <output_subdir>

# 3. Secondary check: numeric comparison of human-readable outputs
python3 compare_results.py --no-intermediates <output_subdir> toml_reference

# 4. Cross-validation: compare toml_reference fitted.toml against v2008 .var
python3 validate_variance.py
```

`clclo2` is excluded by default; pass `--all` to include it.  All four steps
should report 0 failures.

### What each step checks

**Step 2 — `compare_toml_outputs.py` (primary):**
Diffs `mol.fitted.toml` and `mol.catalog.toml` directly against `toml_reference/`.
These files contain the complete numerical output at full double precision:
parameters, full variance matrix, all catalog transitions.  Any change to computed
results is caught here.

**Step 3 — `compare_results.py --no-intermediates` (secondary):**
Numerically compares the human-readable text outputs, which contain content not
present in the TOML files:
- `mol.fit` — per-iteration log (Marquardt parameter, RMS per iteration,
  trust-region changes); `mol.fitted.toml` captures only the final state
- `mol.out` — energy level table, per-transition listing, partition function
- `mol.cat` — catalog lines extracted from `mol.catalog.toml` (redundant with
  step 2, but confirms the extraction is correct)

`--no-intermediates` skips `.par` and `.var`, which TOML mode does not produce.

**Step 4 — `validate_variance.py` (cross-validation):**
Compares `toml_reference/mol.fitted.toml` directly against `v2008_results/mol.var`
for every molecule, checking that the TOML serialization is correct relative to the
independent v2008 baseline:
- Parameter values: relative tolerance ≤ 5×10⁻¹⁴  (`%21.13E` format bound)
- Parameter errors: relative tolerance ≤ 5×10⁻⁷   (`%15.6E` format bound)
- Variance matrix:  absolute tolerance ≤ 5×10⁻⁸    (`%10.7f` format bound)

This is a fixed check on `toml_reference/` rather than `<output_subdir>`, so it
only needs to be re-run after regenerating `toml_reference/`. It's there to ensure the toml_reference
is still correct with respect to the v2008 baseline. Tests against the toml_reference can rely on
`compare_toml_outputs` alone.

### What the TOML path produces

| Program | File | Notes |
|---------|------|-------|
| spfit | `mol.fit`          | Human-readable fit summary; written in both modes |
| spfit | `mol.fitted.toml`  | Full-precision parameters and variance (TOML mode only) |
| spcat | `mol.out`          | Human-readable catalog listing; written in both modes |
| spcat | `mol.catalog.toml` | All catalog transitions (TOML mode only) |

spcat does **not** write a separate `.cat` file in TOML mode — that data lives
inside `mol.catalog.toml`.  `run_toml_tests.py` extracts the `cat_lines` array
and writes it as `mol.cat` for use by `compare_results.py`.  The `.par`, `.var`,
and `.bak` files produced by the legacy path are not written in TOML mode.

### TOML input files

`mol.toml` and `mol.dipoles.toml` are permanent TOML-format inputs stored in each
molecule directory.  On the first run they are generated from the legacy
`.par`/`.lin` and `.int` files via the pickett library and saved for all subsequent
runs.  `run_tests.py` (legacy path) skips `.toml` files when copying to its temp
directory, so the two paths remain independent.

---

## Precision differences: TOML vs legacy

The TOML path stores the variance matrix at full double precision (~17 significant
figures) in `mol.fitted.toml`.  The legacy `.var` format normalizes each row `i`
of the variance by `1/erpar[i]`, writes it as `%10.7f` (7 decimal places), and
`getvar` multiplies back — passing only ~7 significant figures of variance to
`calerr`.

`calerr` computes the ERR column (field 2) of each `.cat` line as a quadratic
form in the variance matrix.  With higher-precision variance, the TOML path
produces slightly more accurate ERR values — differing from v2008 by 1 ULP in
the last printed digit for 17 lines across h2s, o2_1and3, and ch3oh.  Every
other field (frequency, log-intensity, QN columns) is bit-identical to v2008.

`toml_reference/` captures the TOML path's full-precision ERR values and serves
as the regression baseline.  To regenerate it (after a code change that
legitimately alters numerical results):

```sh
uv run --project ../python python3 run_toml_tests.py --all --regenerate-reference toml_reference
python3 compare_results.py --no-intermediates toml_reference v2008_baseline # recheck results based upon legacy files including .cat
python3 validate_variance.py   # re-run cross-validation against v2008_baseline for .par files
```

---

## Notes

- Both test runners copy inputs to a temporary directory before running; originals
  are never overwritten.
- `v2008_results/` is the authoritative reference for the legacy path; any code
  change must preserve it.
- `toml_reference/` is the authoritative reference for the TOML path.
- `reference_outputs/` files were downloaded from the CDMS website and may differ
  slightly from `v2008_results/` (different compiler/era).
- `clclo2` (general_interactions) is excluded from default runs — it is very slow.
