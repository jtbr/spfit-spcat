# pickett Python package

Python bindings for [SPFIT/SPCAT](../README.md) via [nanobind](https://github.com/wjakob/nanobind).

## Prerequisites

- CMake ≥ 3.15  
- A C++17 compiler  
- Python ≥ 3.8 with `pip` / `uv`

## Build and install

```sh
# From the project root, using uv
uv venv && source .venv/bin/activate
uv pip install ".[test]" ./python
```

Or with plain pip:

```sh
pip install "./python[test]"
```

The build compiles the bundled `dblas.c` with `-ffp-contract=off` to preserve
bit-identical numerical results with the v2008 reference baseline.  Do not
override these flags or link an external BLAS.

## Usage

```python
import pickett

# Fit spectroscopic parameters
fit = pickett.fit_files("path/to/molecule")  # reads .par + .lin
print(f"RMS = {fit.xsqbest:.6f}, B = {fit.par[0]:.4f} MHz")

# Generate catalog (no disk output)
cat = pickett.cat_files("path/to/molecule")  # reads .var + .int
for line in cat.cat_lines[:3]:
    print(line)
print(f"Q(300 K) ≈ {dict(zip(cat.temp, cat.qsum)).get(300.0, '?')}")
```

For the low-level path (access to raw session objects and inputs):

```python
session = pickett.FitSession("mol.par", "mol.lin", engine="spinv")
out = session.run()  # may only be called once

session = pickett.CatSession("mol.int", "mol.var", engine="dpi")
out = session.run()
```

## Running tests

```sh
pytest python/tests
```

## Numerical reproducibility

`-ffp-contract=off` disables FMA and `dblas.c` replaces OpenBLAS.  Together
they prevent sign-flipped Cholesky columns and ensure the Python API returns
results byte-identical to the CLI on the same machine.  See `CLAUDE.md` for
the full explanation.
