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
