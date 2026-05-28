"""TOML-path regression tests for SPFIT/SPCAT.

For each molecule in the test suite:
  1. Convert mol.par + mol.lin  → mol.toml         (FitInput TOML, spfit input)
  2. Convert v2008_results/mol.var + mol.int
                                → mol.dipoles.toml  (control + dipoles, spcat input)
  3. Run spfit mol   (TOML auto-detect: reads mol.toml)
     → writes mol.fit + mol.fitted.toml
  4. Run spcat mol   (TOML auto-detect: reads mol.fitted.toml + mol.dipoles.toml)
     → writes mol.out + mol.catalog.toml
  5. Extract cat_lines from mol.catalog.toml → write mol.cat
  6. Move mol.fit, mol.cat, mol.out to <output_subdir>

Compare outputs against the reference baseline using compare_results.py:
    python compare_results.py --no-intermediates <output_subdir> [ref_dir]

Requires the pickett Python package:
    pip install ./python        (from repository root)
    # or: uv sync && uv run ... (from python/)
"""

import os
import subprocess
import shutil
import tempfile
import time
import sys

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomllib  # type: ignore[no-redef]
    except ImportError:
        import tomli as tomllib  # type: ignore[no-redef,import]

try:
    import pickett
    from pickett.toml_io import fit_input_to_dict, cat_input_to_dict, _toml_dumps
    _PICKETT_AVAILABLE = True
except ImportError as _e:
    _PICKETT_AVAILABLE = False
    _PICKETT_IMPORT_ERROR = str(_e)

EXE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
spfit_path = os.path.join(EXE_PATH, "spfit")
spcat_path = os.path.join(EXE_PATH, "spcat")

DEFAULT_REF_DIR = "v2008_results"   # source of mol.var used for dipoles conversion
SPFIT_OUTPUTS   = [".fit"]
SPCAT_OUTPUTS   = [".cat", ".out"]


def _write_toml(path: str, d: dict) -> None:
    with open(path, "wb") as f:
        f.write(_toml_dumps(d).encode())


def run_toml_tests(output_subdir_name: str, suite_base_path: str = ".",
                   skip: list = [], ref_dir_name: str = DEFAULT_REF_DIR) -> None:
    original_cwd = os.getcwd()

    if not _PICKETT_AVAILABLE:
        print(f"Error: pickett package not importable — {_PICKETT_IMPORT_ERROR}")
        print("Install with:  uv pip install -e .  (from python/ directory)")
        return

    if not os.path.exists(spfit_path):
        print(f"Error: spfit not found at {spfit_path}")
        return

    if output_subdir_name == ref_dir_name:
        print(f"Error: output subdir '{output_subdir_name}' is the same as the "
              f"reference dir — refusing to overwrite baseline results.")
        return

    print(f"Starting TOML-path tests for suite in {suite_base_path}")
    print(f"Outputs will be placed in '{output_subdir_name}' subdirectories.")
    print(f"Reference (dipoles source): '{ref_dir_name}'")
    print("-" * 50)

    for category_name in sorted(os.listdir(suite_base_path)):
        category_path = os.path.join(suite_base_path, category_name)
        if not os.path.isdir(category_path):
            continue

        print(f"\nProcessing Category: {category_name}")

        for local_basename in sorted(os.listdir(category_path)):
            mol_path     = os.path.join(category_path, local_basename)
            mol_path_abs = os.path.abspath(mol_path)

            if not os.path.isdir(mol_path):
                continue

            if local_basename in skip:
                print(f"  SKIPPING: {local_basename} (in skip list; use --all to include)")
                continue

            print(f"  Processing Molecule: {local_basename}")

            try:
                with tempfile.TemporaryDirectory() as tmp:
                    # Copy molecule input files into temp dir (writable)
                    for fn in os.listdir(mol_path):
                        src = os.path.join(mol_path, fn)
                        if os.path.isfile(src):
                            dst = os.path.join(tmp, fn)
                            shutil.copy2(src, dst)
                            os.chmod(dst, 0o644)

                    os.chdir(tmp)

                    output_dir = os.path.join(mol_path_abs, output_subdir_name)
                    os.makedirs(output_dir, exist_ok=True)

                    par_path     = os.path.join(tmp, f"{local_basename}.par")
                    lin_path     = os.path.join(tmp, f"{local_basename}.lin")
                    int_path     = os.path.join(tmp, f"{local_basename}.int")
                    ref_var      = os.path.join(mol_path_abs, ref_dir_name, f"{local_basename}.var")
                    toml_in      = os.path.join(tmp, f"{local_basename}.toml")
                    dipoles_toml = os.path.join(tmp, f"{local_basename}.dipoles.toml")
                    fitted_toml  = os.path.join(tmp, f"{local_basename}.fitted.toml")
                    catalog_toml = os.path.join(tmp, f"{local_basename}.catalog.toml")

                    if not os.path.exists(par_path):
                        print(f"    Skipping: {local_basename}.par not found.")
                        continue

                    # ── Step 1: .par + .lin → mol.toml ───────────────────────
                    try:
                        fi = pickett.parse_fit_files(par_path, lin_path)
                        _write_toml(toml_in, fit_input_to_dict(fi))
                        print(f"    Wrote {local_basename}.toml")
                    except Exception as e:
                        print(f"    ERROR writing {local_basename}.toml: {e}")
                        continue

                    # Remove spfit legacy inputs so spfit sees only TOML input.
                    for legacy_ext in (".par", ".lin"):
                        p = os.path.join(tmp, local_basename + legacy_ext)
                        if os.path.exists(p):
                            os.remove(p)

                    # ── Step 2: v2008_results/.var + .int → mol.dipoles.toml ─
                    have_dipoles = False
                    if os.path.exists(int_path):
                        if not os.path.exists(ref_var):
                            print(f"    WARNING: {ref_dir_name}/{local_basename}.var not found; "
                                  f"skipping dipoles.toml (run legacy tests first to generate it).")
                        else:
                            try:
                                ci = pickett.parse_cat_files(ref_var, int_path)
                                _write_toml(dipoles_toml, cat_input_to_dict(ci))
                                print(f"    Wrote {local_basename}.dipoles.toml")
                                have_dipoles = True
                            except Exception as e:
                                print(f"    WARNING: could not write {local_basename}.dipoles.toml: {e}")
                    else:
                        print(f"    Note: {local_basename}.int not found; skipping dipoles.")

                    # Remove spcat legacy inputs so spcat sees only TOML input.
                    for legacy_ext in (".var", ".int"):
                        p = os.path.join(tmp, local_basename + legacy_ext)
                        if os.path.exists(p):
                            os.remove(p)

                    # ── Step 3: Run TOML-mode spfit ───────────────────────────
                    spfit_ok = False
                    spfit_cmd = [spfit_path, local_basename]
                    print(f"    Running SPFIT (TOML): {' '.join(spfit_cmd)}")
                    try:
                        log = os.path.join(output_dir, f"{local_basename}_spfit.log")
                        r = subprocess.run(spfit_cmd, capture_output=True, text=True,
                                           check=False, timeout=600)
                        with open(log, "w") as lf:
                            lf.write("--- STDOUT ---\n"); lf.write(r.stdout)
                            lf.write("\n--- STDERR ---\n"); lf.write(r.stderr)
                        if r.returncode == 0:
                            print("    SPFIT (TOML) succeeded.")
                            spfit_ok = True
                        else:
                            print(f"    SPFIT (TOML) failed (rc={r.returncode}); see {log}")
                    except subprocess.TimeoutExpired:
                        print("    SPFIT (TOML) timed out.")
                    except Exception as e:
                        print(f"    Error running SPFIT: {e}")

                    # ── Step 4: Run TOML-mode spcat ───────────────────────────
                    spcat_ok = False
                    if spfit_ok and os.path.exists(fitted_toml) and have_dipoles:
                        spcat_cmd = [spcat_path, local_basename]
                        print(f"    Running SPCAT (TOML): {' '.join(spcat_cmd)}")
                        try:
                            log = os.path.join(output_dir, f"{local_basename}_spcat.log")
                            r = subprocess.run(spcat_cmd, capture_output=True, text=True,
                                               check=False, timeout=1200)
                            with open(log, "w") as lf:
                                lf.write("--- STDOUT ---\n"); lf.write(r.stdout)
                                lf.write("\n--- STDERR ---\n"); lf.write(r.stderr)
                            if r.returncode == 0:
                                print("    SPCAT (TOML) succeeded.")
                                spcat_ok = True
                            else:
                                print(f"    SPCAT (TOML) failed (rc={r.returncode}); see {log}")
                        except subprocess.TimeoutExpired:
                            print("    SPCAT (TOML) timed out.")
                        except Exception as e:
                            print(f"    Error running SPCAT: {e}")
                    elif spfit_ok and not os.path.exists(fitted_toml):
                        print(f"    Skipping SPCAT: {local_basename}.fitted.toml was not produced.")
                    elif spfit_ok and not have_dipoles:
                        print(f"    Skipping SPCAT: {local_basename}.dipoles.toml not available.")
                    else:
                        print("    Skipping SPCAT (SPFIT did not succeed).")

                    # ── Step 5: Extract catalog lines → mol.cat ───────────────
                    if spcat_ok and os.path.exists(catalog_toml):
                        try:
                            with open(catalog_toml, "rb") as f:
                                cat_d = tomllib.load(f)
                            cat_file = os.path.join(tmp, f"{local_basename}.cat")
                            lines = cat_d.get("cat_lines", [])
                            with open(cat_file, "w") as f:
                                for line in lines:
                                    f.write(line + "\n")
                            print(f"    Extracted {len(lines)} catalog lines → {local_basename}.cat")
                        except Exception as e:
                            print(f"    WARNING: could not extract catalog lines: {e}")

                    # ── Step 6: Move outputs to output subdir ─────────────────
                    to_move = []
                    if spfit_ok:
                        to_move += [f"{local_basename}{ext}" for ext in SPFIT_OUTPUTS]
                    if spcat_ok:
                        to_move += [f"{local_basename}{ext}" for ext in SPCAT_OUTPUTS]
                    for fn in sorted(set(to_move)):
                        src = os.path.join(tmp, fn)
                        if os.path.exists(src):
                            shutil.move(src, os.path.join(output_dir, fn))
                            print(f"      Moved {fn}")
                        else:
                            print(f"      Output {fn} not found (may be normal if spcat not run).")

            except Exception as e:
                print(f"  Outer error for {local_basename}: {e}")
            finally:
                os.chdir(original_cwd)
            time.sleep(0.1)

        print("-" * 50)
    print("\nAll TOML-path tests finished.")
    print(f"\nTo compare against reference outputs, run:")
    print(f"  python compare_results.py --no-intermediates {output_subdir_name} {ref_dir_name}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Run spfit/spcat in TOML mode for all test-suite molecules.")
    ap.add_argument("output_subdir",
                    help="Subdirectory name for outputs within each molecule directory")
    ap.add_argument("--all", dest="include_all", action="store_true",
                    help="Include clclo2 (very slow, excluded by default)")
    ap.add_argument("--ref-dir", default=DEFAULT_REF_DIR, metavar="DIR",
                    help=f"Reference subdirectory supplying mol.var for dipoles "
                         f"conversion (default: {DEFAULT_REF_DIR})")
    args = ap.parse_args()

    skip = [] if args.include_all else ["clclo2"]
    run_toml_tests(args.output_subdir, skip=skip, ref_dir_name=args.ref_dir)
