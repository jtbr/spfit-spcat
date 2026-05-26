# pickett Python package

Python bindings for [SPFIT/SPCAT](../README.md) — the standard spectral fitting
and catalog-generation programs for rotational spectroscopy — built with
[nanobind](https://github.com/wjakob/nanobind).

## Prerequisites

- CMake ≥ 3.15
- A C++17 compiler (e.g. `gcc`)
- Python ≥ 3.8

## Build and install

```sh
# From the python/ directory — recommended for development
uv sync                    # creates .venv and installs with test extras
uv run pytest              # verify the install

# Or with plain pip (from the python/ directory)
pip install -e ".[test]"
pytest
```

> **Never run `cmake` from the repository root** — that overwrites the
> legacy `Makefile`.  The Python build runs CMake in a private temp directory
> automatically; you don't need to invoke it manually.

### Numerical reproducibility

The extension compiles with **`-ffp-contract=off`** and links the bundled
**`dblas.c`** rather than an external BLAS.  These flags are only required if
you need bit-for-bit identical results to the v2008 reference baseline or to
the `spfit`/`spcat` CLI binaries built the same way; for general spectroscopic
use the results are physically equivalent regardless.

The reason they matter for the test suite: FMA and optimised-BLAS routines can
change floating-point accumulation order at the ULP level, which is enough to
flip the sign of a near-zero Householder reflector column in the QR
factorisation inside `lsqfit.c` and thereby negate entire columns of the
Cholesky factor written to `.var` files.  Both outcomes are correct — the
covariance matrix is the same — but the text file differs.  The
`pyproject.toml` / `CMakeLists.txt` enforce both settings by default.

## Usage

Four workflows are available; all are numerically identical.

### 1. File-free (recommended for new code)

Build the input structs directly — no `.par`/`.lin`/`.var`/`.int` files needed.

```python
import pickett
from pickett import (
    FitInput, Parameter, LineRecord,
    EngineOptions, SpinvOptions, VibState,
    FitSession,
)

fi = FitInput(
    title="CO v=0",
    n_iterations=5,
    parameters=[
        Parameter(id=100, value=57635.968, a_priori_error=1e37, label="B"),
        Parameter(id=200, value=-0.184,    a_priori_error=1e37, label="-D"),
        Parameter(id=300, value=1.7e-9,    a_priori_error=1e37, label="H"),
        Parameter(id=400, value=0.0,       fixed=True,           label="L"),
    ],
    lines=[
        LineRecord(qn=[1, 0], nqn=1, freq=115271.2018, err=0.05),
        LineRecord(qn=[2, 1], nqn=1, freq=230538.0000, err=0.05),
        LineRecord(qn=[3, 2], nqn=1, freq=345795.9899, err=0.05),
        LineRecord(qn=[4, 3], nqn=1, freq=461040.7681, err=0.05),
    ],
    engine_options=EngineOptions(spinv=SpinvOptions(vibs=[
        VibState(knmin=0, knmax=0, symmetric_rotor_quanta=True),
    ])),
)

out = FitSession.from_input(fi).run()
print(f"RMS = {out.xsqbest:.4f}   iterations = {out.itr}")
for v, e in zip(out.par, out.erpar):
    print(f"  {v:21.13E}  ±{e:12.5E}")
```

### 2. Parse legacy files → modify → run (bridge workflow)

Reads existing `.par`/`.lin`/`.var`/`.int` files into the typed structs, where
you can inspect or modify any field before running.

```python
from pickett import parse_fit_files, parse_cat_files, FitSession, CatSession

# Fit: load, freeze one parameter, re-run
fi = parse_fit_files("co_4.par", "co_4.lin")
fi.parameters[2].fixed = True               # freeze H
out = FitSession.from_input(fi).run()

# Catalog: load, cap frequency, run
ci = parse_cat_files("co_4.var", "co_4.int")
ci.control.fqmax = 500.0                    # 500 GHz cap
out = CatSession.from_input(ci).run()
```

### 3. One-shot convenience wrappers (equivalent to the CLI tools)

```python
import pickett

fit_out = pickett.fit_files("path/to/molecule")           # reads .par + .lin
cat_out = pickett.cat_files("path/to/molecule")           # reads .var + .int
cat_out = pickett.cat_files("path/to/var_base",
                            int_path="path/to/int_base")  # separate bases

print(f"RMS = {fit_out.xsqbest:.4f},  B = {fit_out.par[0]:.4f} MHz")
print(f"Q(300 K) = {dict(zip(cat_out.temp, cat_out.qsum)).get(300.0):.4f}")
for line in cat_out.cat_lines[:3]:
    print(line)
```

### 4. TOML file format

Human-readable TOML files as an alternative to the legacy fixed-width ASCII formats.  `mol.toml` plays the role of `mol.par` + `mol.lin`;
`mol.var.toml` + `mol.int.toml` play the role of `mol.var` + `mol.int`.

```python
from pickett import (
    load_fit_input, save_fit_output,
    load_cat_input, save_cat_output,
    FitSession, CatSession,
)

# Fit from TOML input
fi  = load_fit_input("co.toml")
out = FitSession.from_input(fi).run()
save_fit_output(out, fi, "co.var.toml")      # fitted params + variance

# Catalog from TOML input
ci      = load_cat_input("co.var.toml", "co.int.toml")
cat_out = CatSession.from_input(ci).run()
save_cat_output(cat_out, "co.cat.toml")
```

The CLI tools auto-detect TOML mode: if `mol.toml` is present, `spfit mol`
uses it and writes `mol.var.toml`; `spcat mol` uses `mol.var.toml` +
`mol.int.toml` and writes `mol.cat.toml`.

To migrate an existing molecule from legacy files, pass `--toml-out` — the
usual legacy output is written plus a `.var.toml` / `.cat.toml` alongside:

```sh
spfit --toml-out mol
spcat --toml-out mol
```

Lower-level helpers for dict-level serialization are also exported:
`fit_input_to_dict`, `fit_input_from_dict`, `cat_input_to_dict`,
`fit_output_to_dict`, `cat_output_to_dict`.

## Output fields

**`CalFitOutput`** — returned by `FitSession.run()` and `fit_files()`:

| Field      | Description |
|------------|-------------|
| `par`      | Fitted parameter values (same order as input) |
| `erpar`    | Estimated 1σ errors |
| `xsqbest`  | Best RMS of (obs − calc) / err |
| `itr`      | Number of iterations performed |
| `variance` | Packed upper-triangular variance matrix (`nfit*(nfit+1)/2` floats); pass to `save_fit_output` so spcat gets accurate ERR values |

**`CalCatOutput`** — returned by `CatSession.run()` and `cat_files()`:

| Field       | Description |
|-------------|-------------|
| `cat_lines` | JPL `.cat`-format strings, sorted by frequency |
| `egy_lines` | `.egy`-format energy lines (when `iflg` enables `EGYFLG`) |
| `str_lines` | `.str`-format strength lines (when `iflg` enables `STRFLG`) |
| `nline`     | Total catalog lines |
| `temp`      | Partition-function temperature grid (K) |
| `qsum`      | Partition function Q(T) at each `temp` |

## Running tests

Tests require the full repository (`spfit_spcat_test_suite/` lives two levels
up from this directory).  `conftest.py` sets the working directory to the
repository root automatically.

```sh
uv run pytest           # from python/
# or
pytest python/tests     # from the repository root
```

## Full API reference

See [`API.md`](../API.md) for the complete struct-by-struct reference,
including all `FitInput`/`CatInput`/`EngineOptions`/`VibState`/`CatControl`
fields, the `EngineKind.Dpi` path, exception hierarchy, and the C++ API.

For the legacy `.par`/`.var`/`.lin`/`.int`/`.cat` file formats and the
underlying spectroscopic conventions (parameter IDs, dipole IDs, QN encoding),
see [`spinv.md`](../spinv.md) (SPFIT/SPCAT) and [`dpi.md`](../dpi.md) (DPFIT/DPCAT).
