# TOML File Format Reference

`spfit` and `spcat` can now read and write TOML files as an alternative to the legacy
fixed-width ASCII formats (`.par`, `.lin`, `.var`, `.int`).  TOML files are
human-readable, comment-friendly, and structurally explicit — preferred for new
work and for archiving.

For the Python/C++ API that reads and writes these files, see [`API.md`](API.md).
For the underlying spectroscopic conventions and parameter/dipole ID schemes,
see [`spinv.md`](spinv.md) (SPFIT/SPCAT) or [`dpi.md`](dpi.md) (DPFIT/DPCAT).

---

## Contents

- [TOML syntax in brief](#toml-syntax-in-brief)
- [File roles](#file-roles)
- [mol.toml — fit input](#moltoml--fit-input)
- [mol.fitted.toml — fit output / cat input](#molfittedtoml--fit-output--cat-input)
- [mol.dipoles.toml — catalog control and dipoles](#moldiptoml--catalog-control-and-dipoles)
- [mol.catalog.toml — catalog output](#molcatalogtoml--catalog-output)
- [CLI auto-detection](#cli-auto-detection)
- [Round-trip workflow](#round-trip-workflow)

---

## TOML syntax in brief

[TOML](https://toml.io) (Tom's Obvious, Minimal Language) is a general-purpose configuration
file format designed to be easy to read and write.  Full specification:
**<https://toml.io/en/v1.0.0>**.

The constructs used in these files:

```toml
# Plain key = value (string, integer, float, bool, inline array)
title         = "CO v=0"
n_iterations  = 4
xsqbest       = 0.5119
extended_qn   = false
spin_degeneracies = [2, 2, 2, 2]

# [table] — a section of named key/value pairs.
# Dot notation creates nested tables.
[engine_options]
kind = "spinv"

  [engine_options.spinv]
  oblate = false
  # Keys inside this block belong to engine_options.spinv

# [[array of tables]] — each [[...]] block appends one element to an array.
# Used for parameters, lines, vibs, dipoles, etc.

[[parameters]]
id = 100; value = 57635.96804; a_priori_error = 1e35; label = "B"

[[parameters]]
id = 200; value = -0.184; a_priori_error = 1e35; label = "-D"
```

Key points:
- `[table]` defines **one** named section.
- `[[array]]` defines **one element** in a sequence; repeat the header for each
  additional element.
- Nesting: `[a.b]` is the `b` subtable inside `a`; `[[a.b]]` appends to the
  array `a.b`.
- Semicolons (`;`) separate key-value pairs on one line — a TOML extension used
  here for compact parameter listings (the syntax highlighter on this page may not understand this)
- Comments start with `#`

---

## File roles

| File                   | Struct / content             | Written by                   | Read by                        |
|------------------------|------------------------------|------------------------------|--------------------------------|
| `mol.toml`             | `FitInput` (full)            | hand-authored or `parse_fit_files` | `spfit`, `load_fit_input` |
| `mol.fitted.toml`      | Fitted params + variance     | `spfit`, `save_fit_output`   | `spcat`, `load_cat_input`      |
| `mol.dipoles.toml`     | Control settings + dipoles   | hand-authored or `parse_cat_files` | `spcat`, `load_cat_input` |
| `mol.catalog.toml`     | Catalog lines + Q(T)         | `spcat`, `save_cat_output`   | consumer code                  |

`spfit` reads either `mol.toml` (TOML mode) or `mol.par` + `mol.lin` (legacy
mode), but never both; see [CLI auto-detection](#cli-auto-detection).  `spcat`
pairs `mol.fitted.toml` with `mol.dipoles.toml` in the same way it pairs
`mol.var` with `mol.int`.

### Files produced in each mode

In TOML mode, `spfit` writes `mol.fitted.toml` and `mol.fit`; the legacy
`.par`, `.var`, and `.bak` files are **not** written.

In TOML mode, `spcat` writes `mol.catalog.toml` and `mol.out`; the separate
`.cat` file is **not** written — catalog lines are available as the
`cat_lines` array inside `mol.catalog.toml`.  `.egy` and `.str` are still
written when `CatControl.iflg` flags request energy-level or per-component
strength output.

`mol.fit` (fit diagnostics: residuals, rms, convergence) and `mol.out`
(preamble: title, settings, parameter table; then partition function table)
are human-readable reports with no TOML equivalent; they are produced in both
modes.

---

## mol.toml — fit input

Corresponds to `FitInput`; read by `spfit` and `load_fit_input`.

```toml
# ── Global fit settings ──────────────────────────────────────────────────────
title        = "CO v=0 (four lines)"
n_iterations = 4           # max fit iterations
# Optional overrides (defaults shown):
# marquardt_param  = 0.0   # Levenberg-Marquardt initial λ
# max_obs_calc_err = 1e6   # reject |obs-calc|/err above this
# param_err_scale  = 1.0   # scale all a_priori_error values (FRAC in spinv.md)
# freq_scale       = 1.0   # infrared frequency scaling (CAL)
# max_lines        = 32767 # soft cap on input lines
# extended_qn      = false # true → allow up to 10 QN per state (needed for
#                          #   molecules with many coupled spins; set automatically
#                          #   when converting from a legacy .par with negative NLINE)
# nxpar            = 0     # trailing params excluded from special-line fits

# ── Engine options ────────────────────────────────────────────────────────────
[engine_options]
kind = "spinv"   # "spinv" (default) or "dpi"

  [engine_options.spinv]
  # Optional global flags (defaults are 0 / false / ""):
  # oblate          = false  # false = prolate (z=a); true = oblate (z=c)
  # inclusion_flags = 0      # IXX: binary flags for inter-state couplings
  # diag_order      = 0      # 0=energy, 1=full proj., 2=energy/Wang, 3=τ, 4=⟨Kz²⟩
  # phase_flags     = 0      # XOPT: PHASE + 10·NEWLZ + 20·NOFC + 40·G12
  # nam_file        = ""     # parameter-label file override (e.g. "spinl.nam")

    [[engine_options.spinv.vibs]]
    # First VibState — carries all global options above.
    # Required fields for non-trivial molecules:
    knmin  = 0    # minimum K (set both to 0 for linear molecules)
    knmax  = 0    # maximum K
    symmetric_rotor_quanta = true   # true for linear / open-shell (K basis)
    spin_degeneracies = [1]         # 2I+1 or 2S+1 for each spin species
    # Optional per-state fields (defaults shown):
    # stat_weight_axis = 1   # IAX: 1=a, 2=b, 3=c, 4=A/B C₂, 6=3-fold, …
    # iwtpl = 1              # statistical weight for even-parity states
    # iwtmn = 1              # statistical weight for odd-parity states
    # vsym  = 0.0            # vibronic symmetry digits (VSYM; last entry only)
    # esym_weight = 99       # EWT; 99 = ignored

    # [[engine_options.spinv.vibs]]  ← add a second block for a second vib state

# ── Parameters ────────────────────────────────────────────────────────────────
# Each [[parameters]] block is one parameter.
# Decimal IDs follow the schemes in spinv.md.
[[parameters]]
id = 100; value = 57635.96804; a_priori_error = 1e35; label = "B"

[[parameters]]
id = 200; value = -0.184; a_priori_error = 1e35; label = "-D"

[[parameters]]
id = 300; value = 1.7e-9; a_priori_error = 1e35; label = "H"

[[parameters]]
id = 400; value = 0.0; fixed = true; label = "L"
# fixed = true marks a parameter as held (NEGBCD flag).  It must immediately
# follow the "master" parameter it is constrained relative to (see below).

# ── Observed transitions ──────────────────────────────────────────────────────
# Use [[raw_lines]] for lines already in legacy .lin text format (fastest,
# preserves exact formatting for bit-identical round-trips).
# Use [[lines]] for structured LineRecord format.

[[raw_lines]]
" 1 0  115271.2018  5.0000e-02  1.0000"

[[raw_lines]]
" 2 1  230538.0000  5.0000e-02  1.0000"
```

### Parameter fields

| Field            | Default  | Notes |
|------------------|----------|-------|
| `id`             | `0`      | Decimal parameter ID from spinv.md (e.g. 100 = *B*) |
| `value`          | `0.0`    | Absolute value in MHz (or cm⁻¹ if `param_err_scale` triggers wavenumber mode) |
| `a_priori_error` | `1e37`   | 1σ prior; `1e37` = essentially free; ≤ 1e-37 effectively fixes the parameter |
| `fixed`          | `false`  | Hold parameter out of the fit; equivalent to `a_priori_error ≤ 1e-37` but more explicit |
| `label`          | `""`     | Display label; truncated to 10 chars in text output |

### Dependent parameters

A *dependent* parameter is one whose value is constrained to be a fixed
multiple of the preceding independent parameter (the negative-IDPAR feature in
spinv.md, §"Dependent Parameters").  To express this:

1. Set `fixed = true` on the dependent parameter.
2. Place it **immediately after** its master in the `[[parameters]]` list.
3. Give it an **absolute value** — the ratio to the master is computed
   internally.

```toml
[[parameters]]
id = 110010000; value = 59501.34; a_priori_error = 1e35; label = "lambda"

[[parameters]]
id = 110010000; value = 178504.03; fixed = true; label = "3lambda"
# 178504.03 / 59501.34 ≈ 3.0 — stored and retrieved as absolute values
```

Both `mol.toml` and `mol.fitted.toml` store absolute values for dependent
parameters; the 3:1 ratio is implicit, not explicit.

### Lines: raw vs structured

`[[raw_lines]]` stores each line as a string in the legacy `.lin` format — the
simplest path when converting from existing files via `parse_fit_files`.

`[[lines]]` stores structured `LineRecord` objects and is more readable for
hand-authored files:

```toml
[[lines]]
qn = [1, 0]; nqn = 1; freq = 115271.2018; err = 0.05; weight = 1.0

[[lines]]
qn = [2, 1]; nqn = 1; freq = 230538.000; err = 0.05
```

`qn` lists upper-state quanta first, then lower-state: for `nqn=3`, that is
`[J', Ka', Kc', J'', Ka'', Kc'']`.  See [`API.md §LineRecord`](API.md#linerecord)
for the full field description.

---

## mol.fitted.toml — fit output / cat input

Written automatically by `spfit` (TOML mode) or `save_fit_output`.  It merges
the fitted parameters and variance from `CalFitOutput` with the engine options
and labels from the original `FitInput`, producing everything `spcat` needs.

**Do not hand-edit** the fitted parameter values or variance array — rounding
errors will corrupt the line-strength uncertainties.  Labels, `itag`, and
dipole entries live in `mol.dipoles.toml` (see below) and are safe to edit.

```toml
title   = "CO v=0"
itr     = 4           # iterations performed
xsqbest = 0.5119      # RMS (obs-calc)/err at convergence
# extended_qn = true  # present only when the molecule needs > 6 QN per state

[engine_options]
# ... same structure as mol.toml; copied verbatim from the input ...

[[parameters]]
id = 100; value = 57635.968040; error = 3.41e-05; label = "B"
# 'error' is the fitted 1σ; loaded as a_priori_error when this file is
# re-read as a FitInput for a subsequent fit.

[[parameters]]
id = 200; value = -0.184001; error = 5.62e-07; label = "-D"

# Dependent parameters store ABSOLUTE values (same convention as mol.toml):
[[parameters]]
id = 110010000; value = 59501.34; error = 0.012; label = "lambda"
[[parameters]]
id = 110010000; value = 178504.03; fixed = true; error = 0.037; label = "3lambda"
# error for dependent param = |ratio| × error of master

# Packed upper-triangular Cholesky factor of the parameter covariance;
# nfit*(nfit+1)/2 floats, column-major.  Required for accurate ERR values
# in the catalog.
variance = [2.37e-05, 1.12e-05, 4.77e-06, ...]
```

The `variance` field is the same data written to the legacy `.var` file.  When
`mol.fitted.toml` is loaded as a new `FitInput` for a subsequent fit (via
`load_fit_input` or by placing it as `mol.toml`), the variance is used as a
correlated Bayesian prior rather than the independent diagonal defaults — tighter
and more physically meaningful.

### `extended_qn`

Set automatically to `true` by `spfit` for molecules that require more than
`MAXCAT = 6` quantum numbers per state (e.g. internal-rotation species with
multiple nuclear spins, like CH₃OH).  `spcat` reads it from `mol.fitted.toml`
and adjusts the quantum-number format accordingly.

If you are hand-constructing `mol.fitted.toml` for such a molecule, add:
```toml
extended_qn = true
```
at the top level.

---

## mol.dipoles.toml — catalog control and dipoles

Hand-authored; pairs with `mol.fitted.toml` for `spcat`.  Corresponds to the
legacy `.int` file.

```toml
# ── Catalog title (optional; appears as line 1 of mol.out) ───────────────────
title = "CO, v = 0"   # if absent, spcat uses the Hamiltonian title from mol.fitted.toml

# ── Control line (maps to .int line 2) ───────────────────────────────────────
[control]
itag   = 28503     # species tag (6-digit catalog identifier)
qrot   = 108.865   # rotational partition function at tmq
lblk   = 99        # end F quantum — must be non-zero for output!
fqmax  = 10330.0   # max frequency in GHz
thrsh  = -40.0     # log-intensity lower cutoff (STR0)
thrsh1 = -40.0     # secondary cutoff (STR1)
# Optional (defaults shown):
# iflg  = 0        # FLAGS: IRFLG·1000 + OUTFLG·100 + STRFLG·10 + EGYFLG
# inblk = 0        # start F quantum (FBGN)
# tmq   = 300.0    # temperature (K) for intensities
# maxv  = -1       # max v; -1 = unlimited

# ── Dipole moments ────────────────────────────────────────────────────────────
[[dipoles]]
id = 1; value = 0.1101135   # μ_a in Debye
# starts_new_component = true  ← True for the first dipole of each component
#                                 when iflg enables STRFLG=2

[[dipoles]]
id = 3; value = 0.1234      # e.g. μ_c
```

**`lblk` is the most common mistake** when hand-authoring this file: the default
is 0, which suppresses all output.  Set it to the maximum F quantum of interest
(typically `99` for molecules without hyperfine structure, or a larger value if
needed).

Dipole IDs follow the same decimal encoding as parameter IDs; see spinv.md
§"FORMAT of the .int file" for the full scheme.

---

## mol.catalog.toml — catalog output

Written by `spcat` (TOML mode) or `save_cat_output`.  In TOML mode this is
the primary catalog output — the separate `.cat` file is not written.
`cat_lines` holds the same fixed-width JPL-format lines that would appear in
`.cat`; `temp` / `qsum` are the partition function table.

```toml
nline = 95

temp = [300.0, 225.0, 150.0, 75.0, 37.5, 18.75, 9.375, 5.0, 2.725]
qsum = [108.865, 81.649, 54.432, 27.216, 13.608, 6.804, 3.402, 2.040, 1.449]

cat_lines = [
  "  115271.2021  0.0001 -5.0105 2    0.0000  3  28503 101 1           0",
  "  230538.0040  0.0001 -4.4390 2    3.8453  5  28503 101 2           1",
  # ... one line per transition, sorted by frequency ...
]
```

Each string in `cat_lines` is a raw JPL-catalog-format line.  See spinv.md
§"FORMAT of .cat" for the fixed-width field layout.  Partition function values
in `qsum` correspond element-by-element to `temp`.

---

## CLI auto-detection

The `spfit` and `spcat` executables detect TOML mode automatically at startup —
no flags needed:

```sh
# TOML mode: reads mol.toml, writes mol.fitted.toml + mol.fit
# (no .par / .var / .bak written)
spfit mol

# Legacy mode: reads mol.par + mol.lin, writes mol.var + mol.par + mol.fit
# (TOML mode takes precedence if mol.toml exists)
spfit mol

# TOML mode: reads mol.fitted.toml + mol.dipoles.toml,
#            writes mol.out + mol.catalog.toml
# (no separate .cat written; cat lines are inside mol.catalog.toml)
spcat mol
```

Detection logic for `spfit`: if `mol.toml` exists in the working directory,
TOML mode is used; the `.par` and `.lin` files are ignored (a warning is printed
if they also exist).  For `spcat`: if both `mol.fitted.toml` and
`mol.dipoles.toml` exist, TOML mode is used.

### Migrating from legacy files with `--toml-out`

To convert an existing molecule to TOML format, run the legacy path once with
`--toml-out`:

```sh
# Produces mol.var + mol.par (legacy) AND mol.fitted.toml (new)
spfit --toml-out mol

# Produces mol.cat (legacy) AND mol.catalog.toml (new)
spcat --toml-out mol
```

For a full migration including `mol.toml` and `mol.dipoles.toml` input files, use the
Python helpers:

```python
import pickett
from pickett.toml_io import fit_input_to_dict, cat_input_to_dict, _toml_dumps

def write_toml(path, d):
    with open(path, "wb") as f:
        f.write(_toml_dumps(d).encode())

fi = pickett.parse_fit_files("mol.par", "mol.lin")
write_toml("mol.toml", fit_input_to_dict(fi))

ci = pickett.parse_cat_files("mol.var", "mol.int")
write_toml("mol.dipoles.toml", cat_input_to_dict(ci))
```

`--toml-out` can be combined with `--dpi` or `--spinv` in any order:

```sh
spfit --dpi --toml-out mol
```

---

## Round-trip workflow

Full TOML round-trip from legacy files — no hand-editing of `.par`/`.var`:

```sh
# 1. Convert legacy inputs to TOML
python3 -c "
import pickett
from pickett.toml_io import fit_input_to_dict, cat_input_to_dict, _toml_dumps
def wt(p, d):
    open(p,'wb').write(_toml_dumps(d).encode())
wt('mol.toml',         fit_input_to_dict(pickett.parse_fit_files('mol.par','mol.lin')))
wt('mol.dipoles.toml', cat_input_to_dict(pickett.parse_cat_files('mol.var','mol.int')))
"

# 2. Remove legacy files to avoid confusion
rm mol.par mol.lin mol.var mol.int

# 3. Fit (TOML auto-detected)
spfit mol        # writes mol.fit + mol.fitted.toml  (no .par / .var / .bak)

# 4. Catalog (TOML auto-detected)
spcat mol        # writes mol.out + mol.catalog.toml  (no separate .cat)
```

Or in pure Python:

```python
import pickett
from pickett.toml_io import (
    fit_input_to_dict, cat_input_to_dict, _toml_dumps,
    load_fit_input, load_cat_input, save_fit_output, save_cat_output,
)

# Parse → TOML → fit
fi = pickett.parse_fit_files("mol.par", "mol.lin")
with open("mol.toml", "wb") as f:
    f.write(_toml_dumps(fit_input_to_dict(fi)).encode())

fi2 = load_fit_input("mol.toml")
fit_out = pickett.FitSession.from_input(fi2).run()
save_fit_output(fit_out, fi2, "mol.fitted.toml")

# fitted.toml + dipoles → catalog
ci = load_cat_input("mol.fitted.toml", "mol.dipoles.toml")
cat_out = pickett.CatSession.from_input(ci).run()
save_cat_output(cat_out, "mol.catalog.toml")

print(f"Q(300 K) = {dict(zip(cat_out.temp, cat_out.qsum))[300.0]:.4f}")
print(cat_out.cat_lines[:3])
```
