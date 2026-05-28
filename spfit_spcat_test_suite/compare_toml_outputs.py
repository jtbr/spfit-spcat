"""Diff TOML output files against a reference baseline.

Compares mol.fitted.toml and mol.catalog.toml produced by a TOML-path test run
against the files stored in a reference subdirectory (default: toml_reference).

Usage:
    python3 compare_toml_outputs.py <output_subdir> [reference_subdir]

Exit code is the number of differing file pairs (0 = all pass).
"""

import difflib
import os
import sys

TOML_SUFFIXES = [".fitted.toml", ".catalog.toml"]
MAX_DIFF_LINES = 20


def compare_toml_outputs(output_subdir: str, reference_subdir: str = "toml_reference",
                         suite_base_path: str = ".") -> int:
    n_differ = 0

    for category_name in sorted(os.listdir(suite_base_path)):
        category_path = os.path.join(suite_base_path, category_name)
        if not os.path.isdir(category_path):
            continue

        for mol_basename in sorted(os.listdir(category_path)):
            mol_path = os.path.join(category_path, mol_basename)
            if not os.path.isdir(mol_path):
                continue

            out_dir = os.path.join(mol_path, output_subdir)
            ref_dir = os.path.join(mol_path, reference_subdir)

            for suffix in TOML_SUFFIXES:
                fname = mol_basename + suffix
                out_file = os.path.join(out_dir, fname)
                ref_file = os.path.join(ref_dir, fname)

                if not os.path.exists(ref_file):
                    continue  # not in reference — not expected for this molecule
                if not os.path.exists(out_file):
                    print(f"  MISSING  {mol_basename}/{output_subdir}/{fname}")
                    n_differ += 1
                    continue

                with open(ref_file) as f:
                    ref_lines = f.readlines()
                with open(out_file) as f:
                    out_lines = f.readlines()

                if ref_lines == out_lines:
                    print(f"  OK       {mol_basename}/{fname}")
                else:
                    n_differ += 1
                    diff = list(difflib.unified_diff(
                        ref_lines, out_lines,
                        fromfile=f"{reference_subdir}/{fname}",
                        tofile=f"{output_subdir}/{fname}",
                        n=2))
                    print(f"  DIFFER   {mol_basename}/{fname}  ({len(diff)} diff lines)")
                    for line in diff[:MAX_DIFF_LINES]:
                        print("    " + line, end="")
                    if len(diff) > MAX_DIFF_LINES:
                        print(f"    ... ({len(diff) - MAX_DIFF_LINES} more lines)")

    print()
    if n_differ == 0:
        print("All TOML outputs match.")
    else:
        print(f"{n_differ} file(s) differ or are missing.")
    return n_differ


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(
        description="Diff TOML output files against a reference baseline.")
    ap.add_argument("output_subdir",
                    help="Subdirectory containing outputs to compare")
    ap.add_argument("reference_subdir", nargs="?", default="toml_reference",
                    help="Reference subdirectory (default: toml_reference)")
    args = ap.parse_args()
    sys.exit(compare_toml_outputs(args.output_subdir, args.reference_subdir))
