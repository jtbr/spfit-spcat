# SPFIT/SPCAT Library API

This document describes how to use the SPFIT/SPCAT spectroscopy code as a
library, primarily through the **typed-struct (file-free) Python API**.

For the format of the legacy `.par`, `.var`, `.lin`, `.int`, `.cat`, `.str`, and
`.egy` files — including the underlying spectroscopic conventions and the full
parameter / dipole ID schemes — see [`spinv.md`](spinv.md) (SPFIT/SPCAT) and
[`dpi.md`](dpi.md) (DPFIT/DPCAT).  These remain authoritative for the
*content* of the structs documented here; this file documents the *Python /
C++ surface* that exposes them.

---

## Contents

- [Overview](#overview)
- [Quick start: fitting CO without files](#quick-start-fitting-co-without-files)
- [Input structs](#input-structs)
  - [`FitInput`](#fitinput)
  - [`Parameter`](#parameter)
  - [`LineRecord`](#linerecord)
  - [`EngineOptions` / `SpinvOptions` / `VibState`](#engine-options)
  - [`DpiOptions`](#dpioptions)
  - [`CatInput`](#catinput)
  - [`DipoleMoment`](#dipolemoment)
  - [`CatControl`](#catcontrol)
- [Sessions and execution](#sessions-and-execution)
- [Output structs](#output-structs)
  - [`CalFitOutput`](#calfitoutput)
  - [`CalCatOutput`](#calcatoutput)
- [Differences from the legacy file format](#differences-from-the-legacy-file-format)
- [Working with legacy files](#working-with-legacy-files)
- [Exceptions](#exceptions)
- [C++ API](#c-api)
- [Installation and build](#installation-and-build)

---

## Overview

Two equivalent entry-point styles are available; **new code should use the
typed-struct API** (left column below), which avoids any temporary files and
makes inputs straightforward to construct, validate, modify, and serialize:

| Typed-struct API (recommended)                | Legacy file-based API                         |
| --------------------------------------------- | --------------------------------------------- |
| `FitSession.from_input(fi).run()`             | `FitSession(par_file, lin_file).run()`        |
| `CatSession.from_input(ci).run()`             | `CatSession(int_file, var_file).run()`        |
| `parse_fit_files(par, lin) → FitInput`        | `fit_files(base) → CalFitOutput`              |
| `parse_cat_files(var, int) → CatInput`        | `cat_files(base, int_path=…) → CalCatOutput`  |

The legacy file-based path is exposed for compatibility with existing `.par`/
`.lin`/`.int`/`.var` files, but it ultimately funnels through the same typed
structs as the file-free path.  `parse_fit_files` / `parse_cat_files` are the
bridges: they read legacy files and return a `FitInput` / `CatInput` that you
can inspect, modify, and pass straight to `from_input`.  Numerical results are
bit-identical between the two routes.

The CLI executables `spfit` and `spcat` remain available for shell-driven use;
they read and write the same legacy files.

---

## Quick start: fitting CO without files

CO is a closed-shell linear molecule: rotational constant $B$, centrifugal
distortion $-D$, $H$, $L$.  The following fits four parameters to four
observed transitions without touching the filesystem.

```python
import pickett
from pickett import (
    FitInput, Parameter, LineRecord,
    EngineOptions, EngineKind, SpinvOptions, VibState,
    FitSession,
)

# --- 1. Parameters (decimal IDs from spinv.md; builder packs to BCD internally) ---
params = []
for pid, val, label in [
    (100, 57635.968,    "B"),
    (200, -0.184,       "-D"),
    (300,  1.7e-9,      "H"),
    (400,  0.0,         "L"),
]:
    p = Parameter()
    p.id = pid
    p.value = val
    p.a_priori_error = 1e37   # very large → free to fit
    p.label = label
    params.append(p)
params[3].fixed = True        # hold L at zero (equivalent to a_priori_error <= 1e-37)

# --- 2. Observed lines (J upper, J lower) ---
def line(ju, jl, freq_mhz, err=0.05):
    lr = LineRecord()
    lr.qn = [ju, jl] + [0] * 18        # qn is a fixed-width 2*MAXQN=20 int array
    lr.nqn = 1                          # one QN per state (J only)
    lr.freq = freq_mhz
    lr.err = err
    lr.weight = 1.0
    return lr

lines = [
    line(1, 0, 115271.2018),
    line(2, 1, 230538.0000),
    line(3, 2, 345795.9899),
    line(4, 3, 461040.7681),
]

# --- 3. Engine options: one VibState, K=0 (linear), symmetric-rotor quanta ---
vib = VibState()
vib.knmin = 0
vib.knmax = 0
vib.symmetric_rotor_quanta = True      # linear molecule → J-only QN scheme
vib.spin_degeneracies = []             # no resolved nuclear hyperfine for ¹²C¹⁶O

so = SpinvOptions()
so.vibs = [vib]

eo = EngineOptions()
eo.kind = EngineKind.Spinv
eo.spinv = so

# --- 4. Assemble FitInput and run ---
fi = FitInput()
fi.title = "CO v=0 — file-free fit"
fi.n_iterations = 5
fi.parameters = params
fi.lines = lines
fi.engine_options = eo

out = FitSession.from_input(fi).run()
print(f"RMS (obs-calc)/err = {out.xsqbest:.4f}   iterations = {out.itr}")
for v, e in zip(out.par, out.erpar):
    print(f"  {v:21.13E}  ± {e:12.5E}")
```

The same molecule, started from existing legacy files instead:

```python
from pickett import parse_fit_files, FitSession

fi = parse_fit_files("co_4.par", "co_4.lin")  # any modifications happen here
out = FitSession.from_input(fi).run()         # numerically identical
```

---

## Input structs

All input structs are exposed by `pickett._pickett` and re-exported from
`pickett`.  They all support default construction followed by attribute
assignment; **keyword-argument constructors are not bound** (`Parameter(id=100)`
will raise `TypeError`).

### `FitInput`

Top-level container for an SPFIT run.

| Field               | Type                       | Default       | Notes |
|---------------------|----------------------------|---------------|-------|
| `title`             | `str`                      | `""`          | Free-form title |
| `n_iterations`      | `int`                      | `1`           | Max fit iterations (`NITR` in spinv.md) |
| `marquardt_param`   | `float`                    | `0.0`         | Initial Marquardt-Levenberg parameter (`THRESH`) |
| `max_obs_calc_err`  | `float`                    | `1e6`         | Reject lines with $|\text{obs−calc}|/\text{err}$ above this (`ERRTST`) |
| `param_err_scale`   | `float`                    | `1.0`         | `FRAC`; see spinv.md §"FORMAT of the .par file" |
| `freq_scale`        | `float`                    | `1.0`         | Infrared frequency scaling (`CAL`) |
| `max_lines`         | `int`                      | `32767`       | `|NLINE|`; soft cap on number of input lines |
| `extended_qn`       | `bool`                     | `False`       | True → allow up to 10 QN per state (`NLINE` negative in spinv.md) |
| `nxpar`             | `int`                      | `0`           | Number of trailing parameters excluded from special-line fits |
| `engine_options`    | `EngineOptions`            | spinv default | See [Engine options](#engine-options) |
| `parameters`        | `list[Parameter]`          | `[]`          | At least one required |
| `variance`          | `list[float]`              | `[]`          | Packed upper-triangular Cholesky decomposition of the prior covariance; empty → diagonal default from each parameter's `a_priori_error` |
| `lines`             | `list[LineRecord]`         | `[]`          | Observed transitions |

The `variance` array, when supplied, holds the lower-left triangle of the
prior parameter covariance Cholesky factor, packed by column: for `nfit`
floated parameters, column `j` contributes `j+1` doubles (so total length is
`nfit*(nfit+1)/2`).  This matches the layout written to the `.var` file.

### `Parameter`

One row in the parameter list.

| Field            | Type    | Default     | Notes |
|------------------|---------|-------------|-------|
| `id`             | `int`   | `0`         | Decimal parameter identifier (e.g. `100` = $B$).  The builder packs this into the BCD `IDPAR` format internally — pass the value you would read from spinv.md's tables, **not** a pre-encoded BCD string. |
| `value`          | `float` | `0.0`       | Absolute value, in MHz unless engine options enable wavenumbers |
| `a_priori_error` | `float` | `1e37`      | 1σ prior; ≤ 1e-37 effectively fixes the parameter |
| `fixed`          | `bool`  | `False`     | True → excluded from fit (sets NEGBCD flag on the BCD-encoded ID) |
| `label`          | `str`   | `""`        | Display label; truncated to 10 chars in any text output |

Set `fixed=True` *or* `a_priori_error ≤ 1e-37` to hold a parameter; the two are
equivalent for the fit but `fixed=True` is clearer in code review.

To express a *dependent* parameter (one whose value follows the previous
parameter in fixed ratio — the negative-IDPAR feature from spinv.md), set
`fixed=True` and order it immediately after its parent in `parameters`.

### `LineRecord`

One observed transition.

| Field      | Type        | Default | Notes |
|------------|-------------|---------|-------|
| `qn`       | `list[int]` | 20 × 0  | Quantum numbers, **flat array of exactly 2·MAXQN = 20 ints**: upper-state quanta first, then lower-state.  Assigning a shorter list raises `TypeError`. |
| `nqn`      | `int`       | `0`     | Number of quanta *per state* that are actually meaningful in `qn` (so 1 ≤ nqn ≤ MAXQN; total slots used = 2·nqn) |
| `freq`     | `float`     | `0.0`   | Observed frequency in MHz (or cm⁻¹ if `err < 0`, matching the legacy `.lin` sign convention) |
| `err`      | `float`     | `1e-7`  | Measurement uncertainty (MHz) |
| `weight`   | `float`     | `1.0`   | Line weight within a blend (`WT`; auto-normalised) |
| `blend_tag`| `str`       | `""`    | Non-empty triggers blend grouping for adjacent lines sharing the tag |

The `qn` layout follows the .lin convention: `qn[0..nqn-1]` are the upper-state
quanta in the order dictated by the quantum-number format (see spinv.md
"Format of Quantum Numbers"), and `qn[MAXQN..MAXQN+nqn-1]` are the lower-state
quanta.  In practice you build the full 20-element list each time, e.g.:

```python
lr.qn = [J_upper, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         J_lower, 0, 0, 0, 0, 0, 0, 0, 0, 0]
lr.nqn = 1
```

Or more compactly for the linear/J-only case:

```python
lr.qn = [J_upper, J_lower] + [0] * 18  # nqn=1 → only slot 0 is read for each state
lr.nqn = 1
```

The latter form works because the formatter ignores positions past `nqn` (and
past `MAXQN` for the lower-state half), but is order-dependent for `nqn ≥ 2`.

### Engine options

```
EngineOptions
├── kind: EngineKind            # EngineKind.Spinv (default) | EngineKind.Dpi
├── spinv: SpinvOptions         # used when kind == Spinv
└── dpi: DpiOptions             # used when kind == Dpi
```

#### `SpinvOptions`

Maps to the first SPINV "option line" of a `.par` file (`CHR, SPIND, NVIB,
KNMIN, KNMAX, IXX, IAX, WTPL, WTMN, VSYM, EWT, DIAG, XOPT` in spinv.md), but
with spectroscopist-friendly field names and the multi-line / sentinel logic
handled by the builder.

| Field              | Type             | Default       | spinv.md name   | Notes |
|--------------------|------------------|---------------|-----------------|-------|
| `inclusion_flags`  | `int`            | `0`           | `IXX`           | Binary flags for inter-state couplings: bit 0 = no ΔN ≠ 0, bit 1 = no ΔJ, etc.  Usually 0. |
| `diag_order`       | `int`            | `0`           | `DIAG`          | Eigenvalue ordering: 0=energy, 1=full projection, 2=energy within Wang, 3=τ=Ka−Kc, 4=⟨K_z²⟩, 5=diag |
| `phase_flags`      | `int`            | `0`           | `XOPT`          | Packed: `PHASE + 10·NEWLZ + 20·NOFC + 40·G12` (see spinv.md) |
| `oblate`           | `bool`           | `False`       | sign of `NVIB`  | True → oblate (z=c); False → prolate (z=a) |
| `nam_file`         | `str`            | `""`          | `CHR`           | Override of the parameter-label file (default `sping.nam`) |
| `vibs`             | `list[VibState]` | `[]`          | one card each   | One entry per vibrational/electronic state; **at least one required** |

#### `VibState`

One vibrational/electronic state's portion of an option line.  `vibs[0]` sets
defaults that propagate to later entries (matching legacy "first option line
sets defaults" semantics).

| Field                     | Type        | Default | spinv.md name      | Notes |
|---------------------------|-------------|---------|--------------------|-------|
| `index`                   | `int`       | `0`     | mag. of `NVIB` (later cards) | 0-based vibrational state index this card configures |
| `knmin`                   | `int`       | `0`     | `KNMIN`            | Minimum K |
| `knmax`                   | `int`       | `359`   | `KNMAX`            | Maximum K; **set both to 0 for linear molecules** |
| `stat_weight_axis`        | `int`       | `1`     | `IAX`              | 1=a, 2=b, 3=c, 4=A/B C₂, 6=3-fold, 7/8=C₄, 9=5-fold, 10/11=C₆.  Negative → use $I_\text{tot}$ basis. |
| `iwtpl`                   | `int`       | `1`     | `WTPL`             | Statistical weight for even-parity states |
| `iwtmn`                   | `int`       | `1`     | `WTMN`             | Statistical weight for odd-parity states |
| `vsym`                    | `float`     | `0.0`   | `VSYM`             | Vibronic symmetry digits (see spinv.md).  Only meaningful on the last entry — the builder injects the `-1.0` continuation sentinel automatically. |
| `esym_weight`             | `int`       | `99`    | `EWT`              | E-symmetry weight.  Default 99 = ignored.  Negative on `vibs[0]` triggers `EWTFAC=1000` mode. |
| `spin_degeneracies`       | `list[int]` | `[]`    | `\|SPIND\|` digits | Spin degeneracy ($2I+1$ or $2S+1$) for each spin species, listed in the order they appear in the QN scheme.  Empty → no spin species. |
| `symmetric_rotor_quanta`  | `bool`      | `False` | sign of `SPIND`    | True for linear molecules and open-shell systems (K-quantum-number basis). |

The `SPIND` field of the legacy format encodes both the sign (QN convention)
and the magnitude (multi-digit spin degeneracies); here those are split into
`symmetric_rotor_quanta` and `spin_degeneracies`, and the builder reassembles
them.  Current limit: individual `spin_degeneracies` values must be in [1, 9]
(single-digit "isbig=0" BCD packing); the multi-digit "isbig=1" path is not
yet exposed.

##### Recipes

```python
# Closed-shell linear molecule (CO, HCN):
vib.knmin = vib.knmax = 0
vib.symmetric_rotor_quanta = True
vib.spin_degeneracies = []

# Closed-shell asymmetric top (H₂O, SO₂):
# leave knmin/knmax at their defaults (0, 359)
vib.spin_degeneracies = [1]    # spin-0 placeholder
vib.stat_weight_axis = 1       # weight on a-axis

# Open-shell linear (O₂, NO):
vib.knmin = vib.knmax = 0
vib.symmetric_rotor_quanta = True
vib.spin_degeneracies = [3]    # 2S+1 for triplet (electron spin)

# Symmetric top with C₃ᵥ (CH₃CN-type) and 3 equivalent H nuclei:
vib.stat_weight_axis = 6       # 3-fold
vib.iwtpl, vib.iwtmn = 2, 1    # A vs E weights
vib.spin_degeneracies = [4]    # I_tot=3/2 from three I=½ → 2·(3/2)+1 = 4
```

For two vibrational states sharing global options (e.g. ar-so2), append a
second `VibState` with `index=1`; the first carries the global `SpinvOptions`
fields and the second only its per-state overrides.

### `DpiOptions`

Minimal DPI engine options (used when `EngineOptions.kind == EngineKind.Dpi`).

| Field    | Type | Default | Notes |
|----------|------|---------|-------|
| `isdgn`  | `int`| `1`     | Spin degeneracy |
| `nvib`   | `int`| `1`     | Number of vibrational states |

See [`dpi.md`](dpi.md) for the underlying DPFIT/DPCAT conventions.  The DPI
engine is exercised by `--dpi` from the CLI; the typed-struct API only
exposes the two fields above (no rich per-state structure yet).

### `CatInput`

Top-level container for an SPCAT run.

| Field             | Type                   | Default       | Notes |
|-------------------|------------------------|---------------|-------|
| `title`           | `str`                  | `""`          | |
| `control`         | `CatControl`           | defaults      | See below |
| `dipoles`         | `list[DipoleMoment]`   | `[]`          | At least one usually required |
| `engine_options`  | `EngineOptions`        | spinv default | Same struct as `FitInput` |
| `parameters`      | `list[Parameter]`      | `[]`          | Same as `FitInput` — the fitted Hamiltonian parameters |
| `variance`        | `list[float]`          | `[]`          | Packed upper-triangular variance matrix (same shape as `FitInput.variance`) |

### `DipoleMoment`

Single row in the dipole list.

| Field                     | Type   | Default | Notes |
|---------------------------|--------|---------|-------|
| `id`                      | `int`  | `0`     | Decimal dipole identifier (`IDIP` in spinv.md §"FORMAT of the .int file") |
| `value`                   | `float`| `0.0`   | Dipole moment in Debye (or unit consistent with parameters) |
| `starts_new_component`    | `bool` | `False` | True → applies the `NEGBCD` flag (begins a new component group when `STRFLG=2`) |

### `CatControl`

Maps to .int line 2 (`FLAGS, TAG, QROT, FBGN, FEND, STR0, STR1, FQLIM, TEMP, MAXV`):

| Field    | Type   | Default   | spinv.md name | Notes |
|----------|--------|-----------|---------------|-------|
| `iflg`   | `int`  | `0`       | `FLAGS`       | `IRFLG·1000 + OUTFLG·100 + STRFLG·10 + EGYFLG` |
| `itag`   | `long` | `999`     | `TAG`         | Catalog species tag |
| `qrot`   | `float`| `1000.0`  | `QROT`        | Rotational partition function at `tmq` |
| `inblk`  | `int`  | `0`       | `FBGN`        | Beginning F quantum (rounded up) |
| `lblk`   | `int`  | `0`       | `FEND`        | Ending F quantum |
| `thrsh`  | `float`| `-100.0`  | `STR0`        | Log-strength cutoff |
| `thrsh1` | `float`| `-100.0`  | `STR1`        | Secondary log-strength cutoff |
| `fqmax`  | `float`| `9999.99` | `FQLIM`       | Maximum frequency in GHz |
| `tmq`    | `float`| `300.0`   | `TEMP`        | Temperature for intensities (K) |
| `maxv`   | `int`  | `-1`      | `MAXV`        | Maximum v; -1 = unlimited |

To enable `.str` output (populates `CalCatOutput.str_lines`), add `10` to
`iflg` (i.e. `STRFLG=1`).  To enable `.egy` output, add the `EGYFLG` digit
(1 = energies, 2 = + derivatives, 3 = + eigenvectors, etc.).  To work in
wavenumbers throughout, add `1000` (`IRFLG=1`).

---

## Sessions and execution

```python
class FitSession:
    @staticmethod
    def from_input(fi: FitInput) -> FitSession: ...      # file-free
    def __init__(self, par_file: str, lin_file: str,
                 engine: str = "spinv"): ...             # legacy files
    def run(self) -> CalFitOutput: ...                   # call exactly once

class CatSession:
    @staticmethod
    def from_input(ci: CatInput) -> CatSession: ...
    def __init__(self, int_file: str, var_file: str,
                 engine: str = "spinv"): ...
    def run(self) -> CalCatOutput: ...
```

`run()` may only be called **once per session**.  Build a fresh session for
each run; this matches the C++ requirement that each `CalFit` / `CalCat`
instance own a single-use `CalculationEngine`.

The `engine` argument of the file-based constructors selects between the
SPINV and DPI engines; the file-free path uses `EngineOptions.kind` instead.

---

## Output structs

### `CalFitOutput`

| Field      | Type          | Description |
|------------|---------------|-------------|
| `par`      | `list[float]` | Fitted parameter values, in the same order as `FitInput.parameters` |
| `erpar`    | `list[float]` | Estimated 1σ errors for each parameter |
| `xsqbest`  | `float`       | Best RMS of (obs − calc) / err across all iterations |
| `itr`      | `int`         | Number of iterations actually performed |

The "extra" fields that `CalFitIO::writeOutput` uses to regenerate text files
(`title_for_output`, `var_final_for_output`, etc.) are intentionally not
exposed in Python — they're only useful for round-tripping through the legacy
file format, which the typed-struct API replaces.

### `CalCatOutput`

| Field       | Type          | Description |
|-------------|---------------|-------------|
| `nline`     | `int`         | Total catalog lines generated |
| `cat_lines` | `list[str]`   | `.cat`-format lines, **sorted by frequency** (matches the CLI's `sortn` step).  Each line is the raw fixed-width JPL catalog format described in spinv.md §"FORMAT of .cat" |
| `egy_lines` | `list[str]`   | `.egy`-format energy lines; populated only when `iflg` enables `EGYFLG` |
| `str_lines` | `list[str]`   | `.str`-format strength lines; populated only when `iflg` enables `STRFLG` |
| `ntemp`     | `int`         | Number of points in the partition-function table |
| `temp`      | `list[float]` | Temperatures (K), `ntemp` entries.  By default this is a descending grid starting at 1000 K. |
| `qsum`      | `list[float]` | Partition function Q(T) at each `temp` entry |

`cat_lines` example (CO J=1→0):

```
  115271.2021  0.0001 -5.0105 2    0.0000  3  28503 101 1           0
  └─freq MHz─┘ ├err┘ ├─lgInt┘ │ ├──ELO──┘ │ ├─tag┘ │ │ └upper Q─┘ └lower Q─┘
              (MHz)  (nm²·MHz)│            │       │ └ NQN
                              DR           GUP     QNFMT
```

Field meanings are defined in spinv.md §"FORMAT of .cat".  Quantum numbers
use Pickett's hybrid base-10/letter encoding for magnitudes above 99 (e.g.
`100` → `A0`, `-10` → `a0`); parse `cat_lines` accordingly.

Useful partition-function lookup:

```python
q_at = dict(zip(c.temp, c.qsum))
print(f"Q(300 K) = {q_at[300.0]:.4f}")
```

---

## Differences from the legacy file format

| Legacy file behaviour                                          | Typed-struct API                                |
|----------------------------------------------------------------|-------------------------------------------------|
| `IDPAR` encoded as BCD digits packed in a `bcd_t` array        | `Parameter.id` is a plain `int`; builder packs to BCD |
| `IDIP` same as above for dipoles                               | `DipoleMoment.id` is a plain `int`              |
| `SPIND` field combines sign + multi-digit spin degeneracies    | Split into `VibState.symmetric_rotor_quanta` and `VibState.spin_degeneracies` |
| `NVIB` sign indicates oblate/prolate                           | `SpinvOptions.oblate` (bool)                    |
| `VSYM = -1.0` sentinel signals "more option lines follow"      | Builder injects this automatically; user only sets the final state's actual `vsym` |
| `IXX`, `IAX`, `EWT`, `DIAG`, `XOPT` short field names          | Renamed to `inclusion_flags`, `stat_weight_axis`, `esym_weight`, `diag_order`, `phase_flags` |
| Negative IDPAR flag means "dependent on previous parameter"    | `Parameter.fixed = True` on a parameter immediately after its parent |
| `a_priori_error ≤ 1e-37` fixes a parameter                     | Either `a_priori_error ≤ 1e-37` *or* `Parameter.fixed = True` |
| `.lin` 12-int QN field; sign of err encodes MHz vs cm⁻¹        | `LineRecord.qn` is 2·MAXQN=20 ints; `err < 0` still encodes wavenumber units |
| `FLAGS = IRFLG·1000 + OUTFLG·100 + STRFLG·10 + EGYFLG`         | `CatControl.iflg` — same packed encoding |
| `.cat` text uses Pickett's letter encoding for \|QN\| > 99     | Same encoding; `cat_lines` are returned as raw strings |
| Lines reformatted on parse / round-trip                        | `parse_fit_files` passes the original `.lin` text through unmodified for bit-identical round-trip |
| Engine selection via `--spinv` / `--dpi` CLI flags             | `EngineOptions.kind = EngineKind.Spinv \| EngineKind.Dpi` |

Things that are **not yet exposed** in the typed-struct API and still
require the legacy file path:

- Multi-digit BCD spin degeneracies (legacy "isbig=1" mode; current limit is digits 1–9).
- DPI engine has only `isdgn` / `nvib` per `DpiOptions`; multi-state DPI runs
  must still use `parse_fit_files` on existing files.

---

## Working with legacy files

The full surface for legacy file handling:

```python
import pickett

# One-shot convenience: read files, run, return results
fit_out = pickett.fit_files("path/to/molecule")              # .par + .lin
cat_out = pickett.cat_files("path/to/molecule")              # .var + .int
cat_out = pickett.cat_files("path/to/var_base",
                            int_path="path/to/int_base")     # separate bases

# Session form (single .run() per instance)
session = pickett.FitSession("co_4.par", "co_4.lin", engine="spinv")
fit_out = session.run()

session = pickett.CatSession("co_4.int", "co_4.var", engine="spinv")
cat_out = session.run()

# Parse to struct, modify, re-run (recommended bridge)
fi = pickett.parse_fit_files("co_4.par", "co_4.lin")
fi.parameters[2].fixed = True                  # freeze H
out = pickett.FitSession.from_input(fi).run()

ci = pickett.parse_cat_files("co_4.var", "co_4.int")
ci.control.fqmax = 500.0                       # cap at 500 GHz
out = pickett.CatSession.from_input(ci).run()
```

`parse_fit_files` / `parse_cat_files` produce structs whose `engine_options`
default to `EngineKind.Spinv`; pass `kind=EngineKind.Dpi` to target the DPI
engine instead.  These parsers preserve the original `.lin`/`.var` text
internally, so round-tripped runs match the legacy-file path to the last bit.

---

## Exceptions

The C++ exception hierarchy is exposed in Python; all classes inherit from
`pickett.CalError`, which itself derives from `Exception`.

| Python class              | Raised when                                              |
|---------------------------|----------------------------------------------------------|
| `pickett.CalError`        | Base class — catch-all for library errors                |
| `pickett.IoError`         | File open/read/write failure                             |
| `pickett.InputError`      | Malformed `.par`/`.lin`/`.int`/`.var` content            |
| `pickett.ValidationError` | API misuse (e.g. `run()` called twice, bad struct shape) |
| `pickett.NumericError`    | Diagonalization or allocation failure                    |

```python
try:
    out = pickett.FitSession.from_input(fi).run()
except pickett.ValidationError as e:
    print(f"Input is malformed: {e}")
except pickett.NumericError as e:
    print(f"Fit failed numerically: {e}")
except pickett.CalError as e:
    print(f"Other library error: {e}")
```

---

## C++ API

The same structs and entry points exist in C++.  Headers of interest:

- `src/api/InputSchema.hpp` — `Parameter`, `LineRecord`, `DipoleMoment`,
  `VibState`, `SpinvOptions`, `DpiOptions`, `EngineOptions`, `CatControl`,
  `FitInput`, `CatInput`.
- `src/api/builders.hpp` — `build_fit_input(fi, engine, logger)` and
  `build_cat_input(ci, engine, logger)`, the bridges to `CalFitInput` /
  `CalCatInput`.
- `src/spfit/CalFit.hpp` / `src/spcat/CalCat.hpp` — the run-once classes
  consumed by the bindings.
- `src/engine/SpinvEngine.hpp` / `src/engine/DpiEngine.hpp` — engine
  implementations of the abstract `CalculationEngine`.
- `src/common/CalError.hpp` — exception hierarchy.

Minimal file-free fit in C++ mirrors the Python example: populate a
`FitInput`, call `build_fit_input(fi, *engine, logger)` to get a
`CalFitInput`, then `CalFit::run(input, output)`.  See
`examples/fit_example.cpp` and `examples/cat_example.cpp` (build with
`-DBUILD_EXAMPLES=ON`).

### `OutputSink` (CalCat only)

C++ `CalCat` writes its text outputs through an abstract `OutputSink`
(`src/spcat/OutputSink.hpp`).  Use `FileSink` for direct file output (CLI
path) or `MemorySink` to capture lines in memory — the Python bindings use
`MemorySink` under the hood.  Each `CalCat` constructor takes four sinks:
`luout` (.out / diagnostics), `lucat` (.cat), `luegy` (.egy), `lustr` (.str).
Mix and match as needed.

### Logger

`src/common/Logger.hpp` provides an abstract `Logger` interface for
diagnostic messages.  Pass a custom subclass to `CalFit` or `CalCat` to
filter or redirect log output:

```cpp
struct QuietLogger : public Logger {
    void log(LogLevel level, const std::string &msg) override {
        if (level >= LogLevel::warn) fprintf(stderr, "%s\n", msg.c_str());
    }
};
QuietLogger quiet;
CalFit calFit(engine, stdout, quiet);
```

---

## Installation and build

### Python

```sh
# system install (rebuilds extension)
pip install ./python

# uv-managed (recommended for development)
cd python && uv pip install -e .
uv run pytest                       # tests; conftest.py sets CWD to repo root
```

By default the build links the system BLAS if available.  For bit-identical
reproduction of the v2008 reference outputs, force the bundled BLAS via:

```sh
pip install ./python --config-settings cmake.args="-DUSE_SYSTEM_BLAS=OFF"
```

The Python `CMakeLists.txt` already enforces `-ffp-contract=off`, the other
half of the FMA / accumulation-order requirement (see `CLAUDE.md` for
rationale).

### C++

```sh
# CMake — always build from build/ (running cmake at the repo root would
# overwrite the legacy Makefile)
mkdir -p build && cd build
cmake .. -DUSE_SYSTEM_BLAS=OFF -DBUILD_EXAMPLES=ON
make
```

### Regression gate

```sh
cd spfit_spcat_test_suite
python3 run_tests.py <output_subdir> && python3 compare_results.py <output_subdir>
```

55 test molecules must reproduce the v2008 reference outputs.

---

## Cross-references

- [`spinv.md`](spinv.md) — SPFIT/SPCAT file formats and spectroscopic conventions (authoritative for parameter IDs, dipole IDs, QN encoding, .cat format).
- [`dpi.md`](dpi.md) — DPFIT/DPCAT file formats (internal rotation).
- [`python/README.md`](python/README.md) — Build prerequisites and packaging.
- [`TASKS.md`](TASKS.md) — Modernization roadmap and status.
