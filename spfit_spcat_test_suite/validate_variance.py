"""Cross-validate mol.fitted.toml against v2008_results/mol.var.

Checks both the fitted parameter values/errors and the variance matrix:

  Parameters: value at %21.13E (~14 sig figs), error at %15.6E (~7 sig figs).
  Variance:   V[i,j] / erpar[i] at %10.7f (~7 sig figs).

A failure indicates a bug in save_fit_output_toml's serialization that the
ERR-column or self-comparison tests cannot detect.

Usage:
    python3 validate_variance.py [suite_base_path]

Exit code is the number of molecules with failures (0 = all pass).
"""

import math
import os
import sys

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomllib  # type: ignore[no-redef]
    except ImportError:
        import tomli as tomllib  # type: ignore[no-redef,import]

TOML_REF   = "toml_reference"
LEGACY_REF = "v2008_results"

# Half a unit in the last place for each format.
TOL_VALUE    = 5e-14   # %21.13E — 13 decimal places in mantissa
TOL_ERROR    = 5e-7    # %15.6E  — 6 decimal places in mantissa
TOL_VARIANCE = 5e-8    # %10.7f  — 7 decimal places


def _row_of(k: int) -> int:
    """Row index i for element k of a packed upper-triangular matrix (column-major).

    Element (i, j) with i <= j is at k = j*(j+1)//2 + i.
    """
    j = int((-1.0 + math.sqrt(1.0 + 8.0 * k)) / 2.0)
    return k - j * (j + 1) // 2


def _parse_var(var_path: str):
    """Parse a .var file.

    Returns (params, erpar_nonfixed, normalized_variance) where:
      params[i]          = (id, value, error) for every parameter in order
      erpar_nonfixed[i]  = error for the i-th non-NEGBCD (id >= 0) parameter
      normalized_variance[k] = V[i,j] / erpar[i] (packed upper-triangular)

    NEGBCD parameters have negative ids and are excluded from the variance matrix.
    All positive-id parameters (including those with small errors like 1e-36) appear
    in the variance.  Some .var files have more than one spin/symmetry header line
    before the parameter list; these are detected by their letter-initial first char.
    """
    with open(var_path) as f:
        lines = f.readlines()

    nfit = int(lines[1].split()[0])

    # Skip spin/symmetry header lines (first non-space char is a letter).
    n_spin = 0
    while lines[2 + n_spin].strip()[0].isalpha():
        n_spin += 1
    param_start = 2 + n_spin

    params = []
    erpar_nonfixed = []
    for idx in range(nfit):
        fields = lines[param_start + idx].split()
        pid   = int(fields[0])
        value = float(fields[1])
        error = float(fields[2])
        params.append((pid, value, error))
        if pid >= 0:
            erpar_nonfixed.append(error)

    n_fitted = len(erpar_nonfixed)
    nvars    = n_fitted * (n_fitted + 1) // 2

    # Variance values: %10.7f, 10 chars wide, 8 per line, no delimiter.
    normalized = []
    for line in lines[param_start + nfit:]:
        row = line.rstrip("\n")
        pos = 0
        while pos + 10 <= len(row):
            chunk = row[pos:pos + 10]
            if chunk.strip():
                normalized.append(float(chunk))
            pos += 10

    return params, erpar_nonfixed, normalized[:nvars]


def _load_toml(fitted_toml_path: str):
    """Return (params, erpar_nonfixed, variance) from mol.fitted.toml.

    params[i]         = (id, value, error) for every parameter in order
    erpar_nonfixed[i] = error for the i-th non-fixed parameter
    variance          = packed upper-triangular variance values
    """
    with open(fitted_toml_path, "rb") as f:
        d = tomllib.load(f)

    raw = d.get("parameters", [])
    params = [
        (p["id"], float(p.get("value", 0.0)),
         float(p.get("error", p.get("a_priori_error", 0.0))))
        for p in raw
    ]
    erpar_nonfixed = [
        float(p.get("error", p.get("a_priori_error", 0.0)))
        for p in raw
        if not p.get("fixed", False)
    ]
    variance = [float(v) for v in d.get("variance", [])]
    return params, erpar_nonfixed, variance


def _reldiff(a: float, b: float) -> float:
    denom = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / denom


def validate(suite_base_path: str = ".") -> int:
    n_fail = 0

    for category in sorted(os.listdir(suite_base_path)):
        cat_path = os.path.join(suite_base_path, category)
        if not os.path.isdir(cat_path):
            continue

        for mol in sorted(os.listdir(cat_path)):
            mol_path = os.path.join(cat_path, mol)
            if not os.path.isdir(mol_path):
                continue

            var_path    = os.path.join(mol_path, LEGACY_REF, f"{mol}.var")
            fitted_path = os.path.join(mol_path, TOML_REF,   f"{mol}.fitted.toml")

            if not os.path.exists(var_path) or not os.path.exists(fitted_path):
                continue

            try:
                var_params, var_erpar, var_norm = _parse_var(var_path)
                toml_params, toml_erpar, toml_var = _load_toml(fitted_path)
            except Exception as e:
                print(f"  ERROR  {mol}: could not parse files: {e}")
                n_fail += 1
                continue

            mol_ok = True

            # ── Parameter values and errors ───────────────────────────────────
            if len(var_params) != len(toml_params):
                print(f"  ERROR  {mol}: param count mismatch "
                      f"(.var={len(var_params)}, toml={len(toml_params)})")
                n_fail += 1
                continue

            max_val_diff = max_err_diff = 0.0
            val_fail = err_fail = False
            for i, ((vid, vval, verr), (tid, tval, terr)) in enumerate(
                    zip(var_params, toml_params)):
                # ids match in magnitude (NEGBCD params have negative id in .var,
                # positive id in TOML).
                if abs(vid) != abs(tid):
                    print(f"  ERROR  {mol}: param {i} id mismatch "
                          f"(.var={vid}, toml={tid})")
                    mol_ok = False
                    break

                vd = _reldiff(tval, vval)
                ed = _reldiff(terr, verr)
                if vd > max_val_diff:
                    max_val_diff = vd
                if ed > max_err_diff:
                    max_err_diff = ed
                if vd > TOL_VALUE:
                    val_fail = True
                if ed > TOL_ERROR:
                    err_fail = True

            if not mol_ok:
                n_fail += 1
                continue

            if val_fail or err_fail:
                mol_ok = False
                n_fail += 1
                if val_fail:
                    print(f"  FAIL   {mol}: max param value reldiff = {max_val_diff:.2e} "
                          f"[tolerance {TOL_VALUE:.0e}]")
                if err_fail:
                    print(f"  FAIL   {mol}: max param error reldiff = {max_err_diff:.2e} "
                          f"[tolerance {TOL_ERROR:.0e}]")

            # ── Variance matrix ───────────────────────────────────────────────
            n_fitted = len(toml_erpar)
            nvars    = n_fitted * (n_fitted + 1) // 2

            if len(var_erpar) != n_fitted:
                print(f"  ERROR  {mol}: non-fixed param count mismatch "
                      f"(.var={len(var_erpar)}, toml={n_fitted})")
                n_fail += 1
                continue

            if len(var_norm) != nvars or len(toml_var) != nvars:
                print(f"  ERROR  {mol}: variance size mismatch "
                      f"(.var={len(var_norm)}, toml={len(toml_var)}, expected={nvars})")
                n_fail += 1
                continue

            max_var_diff = 0.0
            var_fail = False
            for k in range(nvars):
                i  = _row_of(k)
                ep = toml_erpar[i]
                if ep == 0.0:
                    continue
                diff = abs(toml_var[k] / ep - var_norm[k])
                if diff > max_var_diff:
                    max_var_diff = diff
                if diff > TOL_VARIANCE:
                    var_fail = True

            if var_fail:
                mol_ok = False
                n_fail += 1
                print(f"  FAIL   {mol}: max |ΔV_norm| = {max_var_diff:.2e} "
                      f"[tolerance {TOL_VARIANCE:.0e}]")

            if mol_ok:
                print(f"  OK     {mol}: "
                      f"val_reldiff≤{max_val_diff:.1e}  "
                      f"err_reldiff≤{max_err_diff:.1e}  "
                      f"|ΔV_norm|≤{max_var_diff:.1e}  "
                      f"(nfit={n_fitted})")

    print()
    if n_fail == 0:
        print("All parameters and variance matrices validated.")
    else:
        print(f"{n_fail} molecule(s) failed validation.")
    return n_fail


if __name__ == "__main__":
    base = sys.argv[1] if len(sys.argv) > 1 else "."
    sys.exit(validate(base))
