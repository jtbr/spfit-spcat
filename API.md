# SPFIT/SPCAT Library API

This document covers programmatic use of SPFIT/SPCAT from C++ and Python.
For file-format details, see [`spinv.md`](spinv.md) (SPFIT/SPCAT) and [`dpi.md`](dpi.md) (DPFIT/DPCAT).

## Overview

The code can be used in two ways:

- **CLI** — invoke the `spfit` and `spcat` executables with a file base name.
  Files are read from and written to disk exactly as in the original Pickett programs.
- **Library** — instantiate `CalFit` or `CalCat` directly in C++ or Python,
  supply input structs, call `run()`, and read results from output structs.
  No disk I/O is required in the library path (except the scratch file used
  by the block cache inside `CalCat`).

The library path is suitable for embedding in larger analysis workflows,
batch processing, or interactive Python sessions.

---

## C++ API

### Build requirements

The build requires a C++17 compiler and CMake ≥ 3.15.  BLAS selection is
configurable:

```sh
mkdir -p build && cd build

# Use system BLAS (OpenBLAS, MKL, Accelerate, …) — faster on large molecules
cmake ..
make

# Use bundled dblas.c — guarantees bit-identical results with the v2008 reference
cmake .. -DUSE_SYSTEM_BLAS=OFF
make
```

Optimized BLAS libraries and FMA instructions change floating-point
accumulation order at the ULP level, which can flip the sign of near-zero
Householder reflectors in the least-squares solver.  The results are
physically equivalent, but not bit-identical to the v2008 baseline.  Use
`-DUSE_SYSTEM_BLAS=OFF` when exact reproduction of the reference outputs
is required (e.g. regression testing).

To build the C++ examples, add `-DBUILD_EXAMPLES=ON`:

```sh
cmake .. -DBUILD_EXAMPLES=ON
make fit_example cat_example
```

### `CalculationEngine`

`src/engine/CalculationEngine.hpp` declares the abstract interface used by both
`CalFit` and `CalCat` to perform Hamiltonian diagonalization.  Two concrete
implementations are provided:

| Class | Header | Notes |
|---|---|---|
| `SpinvEngine` | `src/engine/SpinvEngine.hpp` | Default; handles asymmetric tops, linear molecules, symmetric tops, electron/nuclear spin interactions (SPFIT/SPCAT). |
| `DpiEngine`   | `src/engine/DpiEngine.hpp`   | Internal-rotation Hamiltonian (DPFIT/DPCAT). |

Engines are heap-allocated and transferred to `CalFit`/`CalCat` via
`std::unique_ptr<CalculationEngine>`.  Each engine instance is single-use —
create a new one for every fit or catalog run.

```cpp
#include "engine/SpinvEngine.hpp"
auto engine = std::make_unique<SpinvEngine>();
```

### `CalFit` — parameter fitting

**Header**: `src/spfit/CalFit.hpp`

#### `CalFitInput`

Populated by `CalFitIO::readInput`; all fields below are set by that call.

| Field | Type | Source | Description |
|---|---|---|---|
| `title` | `std::string` | `.par` line 1 | Molecule title |
| `npar` | `int` | `.par` line 2 | Number of parameters |
| `nitr` | `int` | `.par` line 2 | Requested iteration count |
| `marqp0` | `double` | `.par` line 2 | Initial Marquardt parameter |
| `xerrmx` | `double` | `.par` line 2 | Max allowed (obs−calc)/err |
| `par_initial` | `vector<double>` | `getpar` | Parameter values (MHz or cm⁻¹) |
| `erp_initial` | `vector<double>` | `getpar` | A priori parameter errors |
| `idpar_data` | `vector<unsigned char>` | `getpar` | BCD parameter identifiers |
| `nfit` | `int` | `getpar` | Number of floated parameters |
| `lineData_raw` | `vector<string>` | `.lin` | Observed transition lines |
| `var_initial_from_getvar` | `vector<double>` | `.var` (optional) | Initial variance matrix |

See [`spinv.md`](spinv.md) for the meaning of each `.par` field.

#### `CalFitOutput`

| Field | Type | Description |
|---|---|---|
| `par` | `vector<double>` | Fitted parameter values (same units as input) |
| `erpar` | `vector<double>` | Estimated parameter errors (1σ) |
| `xsqbest` | `double` | Best RMS (obs−calc)/err across all iterations |
| `itr` | `int` | Number of iterations actually performed |

The remaining fields (`title_for_output`, `var_final_for_output`, etc.) are
used by `CalFitIO::writeOutput` to regenerate `.par`/`.var`/`.fit` files;
they can be ignored when only the fitted values are needed.

#### `CalFit` class

```cpp
CalFit(std::unique_ptr<CalculationEngine> &engine,
       FILE *fit_log,
       Logger &logger = Logger::defaultLogger());

void run(const CalFitInput &input, CalFitOutput &output);
```

`engine` ownership is transferred to `CalFit` on construction.
`fit_log` receives the iteration-by-iteration `.fit` text; pass `stdout` or a
`tmpfile()` for silent operation.  `run()` throws a `CalError` subclass on
failure.

#### Full C++ example (`examples/fit_example.cpp`)

```cpp
#include "engine/SpinvEngine.hpp"
#include "spfit/CalFit.hpp"
#include "spfit/CalFitIO.hpp"
#include "common/CalError.hpp"

int main(int argc, char *argv[])
{
    const std::string base = argv[1];          // e.g. "path/to/co_4/co_4"
    auto engine = std::make_unique<SpinvEngine>();

    CalFitInput input;
    if (!CalFitIO::readInput(base + ".par", base + ".lin", input, engine, stdout))
        return EXIT_FAILURE;

    CalFit calFit(engine, stdout);
    CalFitOutput output;
    try {
        calFit.run(input, output);
    } catch (const CalError &e) {
        fprintf(stderr, "Fit failed: %s\n", e.what());
        return EXIT_FAILURE;
    }

    printf("RMS: %.6f  iterations: %d\n", output.xsqbest, output.itr);
    for (size_t i = 0; i < output.par.size(); ++i)
        printf("  [%2zu]  %21.13E  +/-  %12.5E\n", i+1, output.par[i], output.erpar[i]);
}
```

Build and run:

```sh
# from project root
cmake -S . -B build -DBUILD_EXAMPLES=ON && cmake --build build --target fit_example
build/fit_example spfit_spcat_test_suite/diatomic_molecules/co_4/co_4
```

---

### `CalCat` — catalog generation

**Header**: `src/spcat/CalCat.hpp`

#### `OutputSink`

`src/spcat/OutputSink.hpp` provides the abstraction used for all CalCat text
output.  Choose `FileSink` for CLI-style file output, or `MemorySink` to
capture output in memory (required for the Python bindings and useful whenever
you want to process catalog lines without writing files).

```cpp
class OutputSink {
public:
    void printf(const char *fmt, ...);
    virtual void puts(const char *s) = 0;
    virtual FILE *file() const = 0;   // legacy FILE* for engine internals
};

class FileSink : public OutputSink {
    // wraps a FILE* — use for the CLI path
    explicit FileSink(FILE *f);
};

class MemorySink : public OutputSink {
    // accumulates text in a std::string
    std::vector<std::string> drain_lines();  // split on '\n', clear buffer
};
```

`CalCat` takes four `OutputSink*` arguments: `luout` (diagnostic / `.out`),
`lucat` (catalog / `.cat`), `luegy` (energies / `.egy`), `lustr` (strengths /
`.str`).  Any combination of `FileSink` and `MemorySink` is valid.

#### `CalCatInput`

Populated by `CalCatIO::readInput`.

| Field | Type | Source | Description |
|---|---|---|---|
| `title` | `std::string` | `.int` line 1 | Molecule title |
| `itag` | `long` | `.int` line 2 | Species tag (CDMS catalog convention) |
| `qrot` | `double` | `.int` line 2 | Rotational partition function (reference) |
| `thrsh` / `thrsh1` | `double` | `.int` line 2 | Intensity thresholds (log₁₀ nm²·MHz) |
| `fqmax` | `double` | `.int` line 2 | Upper frequency limit (MHz) |
| `dip` | `vector<double>` | `.int` body | Dipole moment components (Debye) |
| `idip` | `vector<unsigned char>` | `.int` body | BCD dipole type identifiers |
| `npar` / `nfit` | `int` | `.var` via `getpar` | Parameter counts |
| `par` | `vector<double>` | `.var` | Hamiltonian parameters |
| `var` | `vector<double>` | `.var` | Variance matrix (packed upper triangular) |

See [`spinv.md`](spinv.md) for the `.int` and `.var` formats in full.

#### `CalCatOutput`

| Field | Type | Description |
|---|---|---|
| `nline` | `long` | Total catalog transitions generated |
| `cat_lines` | `vector<string>` | `.cat` format lines (empty if `FileSink` used) |
| `egy_lines` | `vector<string>` | `.egy` format lines |
| `str_lines` | `vector<string>` | `.str` format lines |
| `ntemp` | `int` | Number of temperature points in partition function table |
| `temp[ntemp]` | `double[]` | Temperatures (K) |
| `qsum[ntemp]` | `double[]` | Partition function Q(T) at each temperature |

`sort_cat_lines()` sorts `cat_lines` by frequency (columns 0–12) then quantum
numbers (columns 55–94), matching the CLI's `sortn` post-processing step.
Call it after draining `MemorySink`s.

#### `CalCat` class

```cpp
CalCat(std::unique_ptr<CalculationEngine> &engine,
       OutputSink *luout, OutputSink *lucat, OutputSink *luegy, OutputSink *lustr,
       Logger &logger = Logger::defaultLogger());

void run(const CalCatInput &input, CalCatOutput &output);
```

#### Full C++ example (`examples/cat_example.cpp`)

```cpp
#include "engine/SpinvEngine.hpp"
#include "spcat/CalCat.hpp"
#include "spcat/CalCatIO.hpp"
#include "spcat/OutputSink.hpp"
#include "common/CalError.hpp"

int main(int argc, char *argv[])
{
    // argv[1] = int_base, argv[2] = var_base
    const std::string int_file = std::string(argv[1]) + ".int";
    const std::string var_file = std::string(argv[2]) + ".var";

    FileSink out_sink(stdout);                         // diagnostic → stdout
    MemorySink cat_sink, egy_sink, str_sink;           // capture .cat/.egy/.str

    auto engine = std::make_unique<SpinvEngine>();
    CalCatInput input;
    if (!CalCatIO::readInput(int_file, var_file, input, engine, &out_sink))
        return EXIT_FAILURE;

    CalCat calCat(engine, &out_sink, &cat_sink, &egy_sink, &str_sink);
    CalCatOutput output;
    try {
        calCat.run(input, output);
    } catch (const CalError &e) {
        fprintf(stderr, "Catalog failed: %s\n", e.what());
        return EXIT_FAILURE;
    }

    output.cat_lines = cat_sink.drain_lines();
    output.egy_lines = egy_sink.drain_lines();
    output.str_lines = str_sink.drain_lines();
    output.sort_cat_lines();   // sort by frequency, like the CLI

    printf("%ld lines generated\n", output.nline);
    for (int i = 0; i < std::min(10, (int)output.cat_lines.size()); ++i)
        printf("%s\n", output.cat_lines[i].c_str());
}
```

Build and run (note: `.var` is the fitted output from spfit):

```sh
cmake -S . -B build -DBUILD_EXAMPLES=ON && cmake --build build --target cat_example
build/cat_example spfit_spcat_test_suite/diatomic_molecules/co_4/co_4 \
                  spfit_spcat_test_suite/diatomic_molecules/co_4/v2008_results/co_4
```

---

### Exception hierarchy

All library errors derive from `CalError` (`src/common/CalError.hpp`):

```
CalError  (std::runtime_error)
├── IoError          — file open/read/write failures
├── InputError       — malformed .par/.lin/.int/.var content
├── ValidationError  — API misuse (wrong sizes, null pointers, re-used sessions)
└── NumericError     — diagonalization failure, memory allocation error
```

Each exception carries a `CalErrorCode` enum for programmatic inspection:

```cpp
try {
    calFit.run(input, output);
} catch (const IoError &e) {
    // file problem
} catch (const NumericError &e) {
    // numerical failure
} catch (const CalError &e) {
    // anything else
}
```

---

### `Logger`

`src/common/Logger.hpp` provides an abstract logging interface.

```cpp
class Logger {
public:
    virtual void log(LogLevel level, const std::string &msg) = 0;
    void debug/info/warn/error(const std::string &msg);     // convenience
    void debug/info/warn/error(const char *fmt, A a, ...);  // printf-style
    static Logger &defaultLogger();   // writes info/warn/error to stdout/stderr
};
```

Pass a custom `Logger` subclass to `CalFit` or `CalCat` to redirect or filter
log messages:

```cpp
struct QuietLogger : public Logger {
    void log(LogLevel level, const std::string &msg) override {
        if (level >= LogLevel::warn)
            fprintf(stderr, "[warn] %s\n", msg.c_str());
    }
};
QuietLogger quiet;
CalFit calFit(engine, stdout, quiet);
```

---

## Python API

### Installation

Prerequisites: a C++17 compiler, CMake ≥ 3.15, and `uv` (or `pip`).

```sh
pip install ./python
```

or, to install into a `uv`-managed environment and run tests:

```sh
uv run --directory python pytest tests/
```

The build compiles the nanobind extension `_pickett` against the project's C++
sources with `-ffp-contract=off` enforced in `python/CMakeLists.txt`.  By
default it links the system BLAS if available; pass
`-DUSE_SYSTEM_BLAS=OFF` to the CMake invocation (via
`pip install ./python --config-settings cmake.args="-DUSE_SYSTEM_BLAS=OFF"`)
for bit-identical reproduction of the v2008 reference outputs.

### High-level API

```python
import pickett

# Fit spectroscopic parameters
fit_out = pickett.fit_files("path/to/molecule")
# fit_out.par      — list[float] of fitted parameter values
# fit_out.erpar    — list[float] of 1σ errors
# fit_out.xsqbest  — float, best RMS (obs−calc)/err
# fit_out.itr      — int, iteration count

# Generate catalog (no files written)
cat_out = pickett.cat_files("path/to/molecule")
# cat_out.cat_lines  — list[str], catalog lines sorted by frequency
# cat_out.egy_lines  — list[str], energy level lines
# cat_out.str_lines  — list[str], intensity/strength lines
# cat_out.nline      — int, total transitions generated
# cat_out.temp       — list[float], temperatures (K)
# cat_out.qsum       — list[float], Q(T) at each temperature

# Engine selection: "spinv" (default) or "dpi"
fit_out = pickett.fit_files("path/to/molecule", engine="spinv")
```

`fit_files(base_path)` reads `base_path + ".par"` and `base_path + ".lin"`.
`cat_files(base_path)` reads `base_path + ".var"` and `base_path + ".int"`.
An optional `int_path` keyword overrides the `.int` base when the intensity
file lives elsewhere.

### Low-level API

For finer control, use the session objects directly.  Each session reads input
files at construction time and may call `run()` exactly once.

```python
import pickett

# --- Fitting ---
session = pickett.FitSession(
    par_file="co_4.par",
    lin_file="co_4.lin",
    engine="spinv",      # default
)
fit_out = session.run()
print(fit_out.xsqbest, fit_out.par)

# --- Catalog ---
session = pickett.CatSession(
    int_file="co_4.int",
    var_file="co_4.var",
    engine="spinv",
)
cat_out = session.run()
q300 = dict(zip(cat_out.temp, cat_out.qsum)).get(300.0)
print(f"Q(300 K) = {q300:.4f}")
for line in cat_out.cat_lines[:5]:
    print(line)
```

### File-free (typed-struct) API

The file-free path lets you construct, modify, or generate input programmatically —
no `.par`, `.lin`, or `.int` files needed.

#### Entry points

| Function / method | What it does |
|---|---|
| `parse_fit_files(par_path, lin_path)` | Parse existing files → `FitInput` struct |
| `parse_cat_files(var_path, int_path)` | Parse existing files → `CatInput` struct |
| `FitSession.from_input(fi)` | Build a fit session from a `FitInput` struct |
| `CatSession.from_input(ci)` | Build a cat session from a `CatInput` struct |

`parse_fit_files` / `parse_cat_files` are useful as a **starting point**: parse
existing files, inspect or modify the resulting struct, then hand it to
`from_input`.  The fit-line content is transferred internally without
reformatting, so the numerical results are bit-identical to the file-based path.

#### Struct overview

```
FitInput
├── title: str
├── n_iterations: int          (default 1)
├── marquardt_param: float     (default 0)
├── max_obs_calc_err: float    (default 1e6)
├── parameters: list[Parameter]
│     Parameter.id              — decimal parameter ID (e.g. 100 = B, 200 = D)
│     Parameter.value           — value in MHz (or cm⁻¹ in cm-wave mode)
│     Parameter.a_priori_error  — 1σ prior; very large → float freely; ≤0 → fix
│     Parameter.fixed           — True → excluded from fit
│     Parameter.label           — display label (≤ 10 chars)
├── lines: list[LineRecord]
│     LineRecord.qn             — list of ints: [Jupper, ...], [Jlower, ...] interleaved
│     LineRecord.nqn            — number of QN pairs actually used
│     LineRecord.freq           — observed frequency (MHz)
│     LineRecord.err            — measurement uncertainty (MHz)
│     LineRecord.weight         — line weight (default 1.0)
└── engine_options: EngineOptions
      EngineOptions.kind        — EngineKind.Spinv (default) or EngineKind.Dpi
      EngineOptions.spinv       — SpinvOptions (see below)
```

`CatInput` is analogous, with `dipoles: list[DipoleMoment]` instead of lines and a
`CatControl` block for catalog settings.

#### SpinvOptions and VibState

Engine options describe the molecular symmetry, quantum numbers, and statistical
weights.  Every fit/cat run needs at least one `VibState` entry.

```python
from pickett import SpinvOptions, VibState

# Linear molecule (e.g. CO, HCN) — single vibrational state
so = SpinvOptions()
so.vibs = [VibState()]
vib = so.vibs[0]
vib.knmin = 0                      # K range → 0,0 signals linear geometry
vib.knmax = 0
vib.symmetric_rotor_quanta = True  # J is the only quantum number; use K-basis
vib.spin_degeneracies = []         # 12C16O: both nuclei have spin 0 → no hyperfine
vib.iwtpl = 1                      # statistical weight for even-parity states
vib.iwtmn = 1                      # statistical weight for odd-parity states

# Closed-shell asymmetric top (e.g. H2O, SO2) — no nuclear spins of interest
so = SpinvOptions()
so.vibs = [VibState()]
vib = so.vibs[0]
# knmin/knmax: defaults (0, 359) include all Ka values → fine for most tops
vib.spin_degeneracies = [1]        # placeholder spin-0 species; gives correct weights
vib.iwtpl = 1
vib.iwtmn = 1

# Symmetric top with C3v symmetry (e.g. CH3CN) — A/E nuclear spin statistics
so = SpinvOptions()
so.vibs = [VibState()]
vib = so.vibs[0]
vib.stat_weight_axis = 6           # 3-fold symmetric top
vib.iwtpl = 2                      # A species weight
vib.iwtmn = 1                      # E species weight
vib.spin_degeneracies = [4]        # three equivalent H with I=1/2: 2*(3/2)+1 = 4
```

Key `VibState` fields:

| Field | Default | Meaning |
|---|---|---|
| `knmin`, `knmax` | `0, 359` | K-quantum-number range. Set both to 0 for linear molecules. |
| `stat_weight_axis` | `1` | Rotation axis for statistical weights: 1=a, 2=b, 3=c; 4=A/B C₂; 6=3-fold; 7=C₄; 10=C₆. Negative → I_tot spin basis. |
| `iwtpl`, `iwtmn` | `1, 1` | Statistical weights for +/− parity states. |
| `vsym` | `0.0` | Vibrational symmetry. Only meaningful on the last (or only) VibState. |
| `esym_weight` | `99` | E-symmetry weight. Default 99 = ignored. Negative → EWTFAC=1000 mode. |
| `spin_degeneracies` | `[]` | Spin degeneracy (2I+1 or 2S+1) for each spin species. `[]` = no spin species. |
| `symmetric_rotor_quanta` | `False` | True for linear molecules and open-shell systems using K quantum numbers. |

Key `SpinvOptions` fields:

| Field | Default | Meaning |
|---|---|---|
| `inclusion_flags` | `0` | Which off-diagonal blocks to include (0 = all). Rarely changed. |
| `diag_order` | `0` | Eigenvalue ordering within Wang sub-blocks (0 = energy order). |
| `phase_flags` | `0` | Packed: PHASE + 10·NEWLZ + 20·NOFC + 40·G12. |
| `oblate` | `False` | True for oblate tops (z=c). Default prolate (z=a). |
| `nam_file` | `""` | Parameter-label file (empty = engine default `sping.nam`). |

#### Example: fit CO (file-free)

CO is a closed-shell linear molecule. Parameters use the standard Pickett ID
scheme: B=100, D=200, H=300, L=400 (centrifugal distortion corrections).

```python
import pickett
from pickett import (
    FitInput, Parameter, LineRecord,
    EngineOptions, EngineKind, SpinvOptions, VibState,
    FitSession,
)

# Spectroscopic parameters (starting values in MHz)
params = [
    Parameter(id=100, value=57635.968,   a_priori_error=1e37),  # B
    Parameter(id=200, value=-0.184,      a_priori_error=1e37),  # -D (sign convention)
    Parameter(id=300, value=1.7e-9,      a_priori_error=1e37),  # H
    Parameter(id=400, value=0.0,         a_priori_error=1e37),  # L (fix at zero)
]
params[3].fixed = True   # L → excluded from fit

# Observed transitions: LineRecord.qn = [J_upper, J_lower, ...] (just J for linear)
# frequencies in MHz (illustrative values from CDMS/JPL)
lines = [
    LineRecord(qn=[1, 0], nqn=1, freq=115271.202, err=0.05),
    LineRecord(qn=[2, 1], nqn=1, freq=230538.000, err=0.05),
    LineRecord(qn=[3, 2], nqn=1, freq=345795.990, err=0.05),
    LineRecord(qn=[4, 3], nqn=1, freq=461040.768, err=0.05),
]

# Engine options: linear molecule → K=0, symmetric-rotor quanta, no nuclear spins
vib = VibState()
vib.knmin = 0
vib.knmax = 0
vib.symmetric_rotor_quanta = True
vib.spin_degeneracies = []

so = SpinvOptions()
so.vibs = [vib]

fi = FitInput()
fi.title = "CO v=0 (file-free example)"
fi.n_iterations = 5
fi.parameters = params
fi.lines = lines
fi.engine_options = EngineOptions()
fi.engine_options.kind = EngineKind.Spinv
fi.engine_options.spinv = so

out = FitSession.from_input(fi).run()
print(f"RMS = {out.xsqbest:.4f}  itr = {out.itr}")
for p, e in zip(out.par, out.erpar):
    print(f"  {p:21.13E}  +/-  {e:12.5E}")
```

#### Parse-and-modify pattern

The most common use of the file-free API is to parse existing files, inspect or
modify the result, and then re-run:

```python
import pickett
from pickett import parse_fit_files, FitSession

# Round-trip: files → struct → run (numerically identical to fit_files())
fi = parse_fit_files("co_4.par", "co_4.lin")
out = FitSession.from_input(fi).run()

# Modify: fix a parameter and re-run
for p in fi.parameters:
    if p.id == 200:          # D constant
        p.a_priori_error = 1e-37  # ≤ 1e-37 effectively fixes it
        break
out2 = FitSession.from_input(fi).run()
print(f"3-param fit RMS = {out2.xsqbest:.4f}")
```

### Exception handling

The C++ exception hierarchy is mapped to Python exceptions of the same names,
all inheriting from `pickett.CalError` (which itself inherits from the built-in
`Exception`):

| Python exception | Raised when |
|---|---|
| `pickett.CalError` | Base class; catch-all for library errors |
| `pickett.IoError` | File open/read/write failure |
| `pickett.InputError` | Malformed `.par`/`.lin`/`.int`/`.var` content |
| `pickett.ValidationError` | API misuse (e.g. calling `run()` twice) |
| `pickett.NumericError` | Diagonalization or allocation failure |

```python
import pickett

try:
    fit_out = pickett.fit_files("missing/molecule")
except pickett.IoError as e:
    print(f"Could not open input files: {e}")
except pickett.CalError as e:
    print(f"Library error: {e}")
```

### Worked Python example

```python
import pickett
from pathlib import Path

BASE = "spfit_spcat_test_suite/diatomic_molecules/co_4/co_4"

# --- Fit ---
fit = pickett.fit_files(BASE)
print(f"Fit converged: RMS={fit.xsqbest:.4f}  iterations={fit.itr}")
for i, (p, e) in enumerate(zip(fit.par, fit.erpar), 1):
    print(f"  [{i:2d}]  {p:21.13E}  +/-  {e:12.5E}")

# --- Catalog ---
# Use the pre-computed .var from the test suite reference outputs
VAR_BASE = "spfit_spcat_test_suite/diatomic_molecules/co_4/v2008_results/co_4"
cat = pickett.cat_files(BASE, int_path=BASE)  # .int from BASE, .var inferred as BASE.var
# or equivalently, override both:
# cat = pickett.CatSession(int_file=BASE+".int", var_file=VAR_BASE+".var").run()

q_by_temp = dict(zip(cat.temp, cat.qsum))
print(f"\nQ(300 K) = {q_by_temp.get(300.0, float('nan')):.4f}")
print(f"Catalog lines: {cat.nline}")
print("First 5 lines:")
for line in cat.cat_lines[:5]:
    print(" ", line)
```

---

## Cross-references

- [`spinv.md`](spinv.md) — File formats for SPFIT/SPCAT: `.par`, `.lin`, `.var`, `.int`, `.cat`, `.egy`, `.str`.
- [`dpi.md`](dpi.md) — File formats for DPFIT/DPCAT (internal rotation).
- [`python/README.md`](python/README.md) — Build prerequisites and packaging details for the Python module.
- [`TASKS.md`](TASKS.md) — Modernization roadmap and status.
