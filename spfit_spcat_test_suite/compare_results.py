import os
import re
import math
import sys

# Tolerances for float comparisons
DEFAULT_REL_TOL = 1e-5  # Relative tolerance
DEFAULT_ABS_TOL = 1e-8  # Absolute tolerance for numbers close to zero

# Regex to find numbers (integers, floats, scientific notation)
# Handles: 123, 123.45, .45, 123., 1e5, -1.2e-3, 1.234E+002
NUMBER_REGEX_STR = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
NUMBER_RE = re.compile(NUMBER_REGEX_STR)

# Regex to find and remove common date/time stamps
# e.g., "Thu Feb 23 16:19:27 2006" or "Mon Jan 01 12:00:00 2000"
DATE_REGEX_STR = r"(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun)\s+(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\s+\d{1,2}\s+\d{2}:\d{2}:\d{2}(?:\.\d+)?\s+\d{4}\b"
DATE_RE = re.compile(DATE_REGEX_STR)

def normalize_line_text(line_text):
    line_no_dates = DATE_RE.sub("[DATE_REMOVED]", line_text)
    line_normalized_equals = re.sub(r'\s*=\s*', '=', line_no_dates)
    return line_normalized_equals

def compare_lines(line1_orig, line2_orig, line_num_str, file_type, rel_tol, abs_tol):
    """
    Compares two lines, ignoring dates and comparing numbers with tolerance.
    Returns True if lines match, False otherwise, along with a description of the difference.
    """
    line1_processed = normalize_line_text(line1_orig)
    line2_processed = normalize_line_text(line2_orig)

    # Handle cases where one line becomes empty after date removal, but other isn't
    if not line1_processed.strip() and not line2_processed.strip():
        return True, ""
    if line1_processed.strip() and not line2_processed.strip():
        return False, f"{line_num_str}: Content mismatch (Ref has content, New is blank after normalization).\n  Ref: '{line1_orig.strip()}'\n  New: '{line2_orig.strip()}'"
    if not line1_processed.strip() and line2_processed.strip():
        return False, f"{line_num_str}: Content mismatch (Ref is blank, New has content after normalization).\n  Ref: '{line1_orig.strip()}'\n  New: '{line2_orig.strip()}'"

    nums1_str = NUMBER_RE.findall(line1_processed)
    nums2_str = NUMBER_RE.findall(line2_processed)

    try:
        nums1 = [float(n) for n in nums1_str]
        nums2 = [float(n) for n in nums2_str]
    except ValueError as e:
        # This might happen if NUMBER_REGEX incorrectly captures something non-numeric
        # or if a non-numeric token that looks like a number (e.g. version "1.0a") appears
        return False, f"{line_num_str}: Error converting extracted 'numbers' to float ({e})."\
              f"\n  Ref line (norm): '{line1_processed.strip()}' (extr: {nums1_str})"\
              f"\n  New line (norm): '{line2_processed.strip()}' (extr: {nums2_str})"

    text_skeleton1 = ' '.join(NUMBER_RE.sub(" _NUM_ ", line1_processed).strip().split())
    text_skeleton2 = ' '.join(NUMBER_RE.sub(" _NUM_ ", line2_processed).strip().split())

    if text_skeleton1 != text_skeleton2:
        return False, f"{line_num_str}: Text structure mismatch."\
                      f"\n          Ref skel: '{text_skeleton1}' (from: '{line1_orig.strip()}')"\
                      f"\n          New skel: '{text_skeleton2}' (from: '{line2_orig.strip()}')"

    if len(nums1) != len(nums2):
        return False, f"{line_num_str}: Different number of numerical values."\
                      f"\n          Ref ({len(nums1)}): {nums1_str} (from: '{line1_orig.strip()}')"\
                      f"\n          New ({len(nums2)}): {nums2_str} (from: '{line2_orig.strip()}')"

    for i, (n1, n2) in enumerate(zip(nums1, nums2)):
        if not compare_floats(n1, n2, rel_tol, abs_tol):
            return False, f"{line_num_str}: Numerical mismatch at value #{i+1}."\
                          f"\n          Ref val: {n1} (from token: {nums1_str[i]})"\
                          f"\n          New val: {n2} (from token: {nums2_str[i]})"\
                          f"\n          Ref line: '{line1_orig.strip()}'"\
                          f"\n          New line: '{line2_orig.strip()}'"
    return True, ""

def compare_floats(f1, f2, rel_tol, abs_tol):
    """Compares two floats within given tolerances."""
    if math.isinf(f1) and math.isinf(f2) and (f1 > 0) == (f2 > 0): return True
    if math.isnan(f1) and math.isnan(f2): return True
    return math.isclose(f1, f2, rel_tol=rel_tol, abs_tol=abs_tol)

def compare_file_pair_final(ref_file_path, new_file_path, file_type, rel_tol, abs_tol):
    """Compares two files, line by line."""
    print(f"    Comparing {file_type} files:\n      Ref: {ref_file_path}\n      New: {new_file_path}")

    if not os.path.exists(ref_file_path):
        print(f"    RESULT: Reference file {os.path.basename(ref_file_path)} not found. Skipping.\n")
        return "ref_missing", []
    if not os.path.exists(new_file_path):
        print(f"    RESULT: New output file {os.path.basename(new_file_path)} not found. Cannot compare.\n")
        return "new_missing", []

    differences = []

    try:
        with open(ref_file_path, 'r', encoding='latin-1') as f_ref, \
             open(new_file_path, 'r', encoding='latin-1') as f_new:

            ref_lines_orig = f_ref.readlines()
            new_lines_orig = f_new.readlines()

            i, j = 0, 0 # Pointers for ref_lines and new_lines

            while i < len(ref_lines_orig) and j < len(new_lines_orig):
                line_num_str = f"Ref L{i+1}/New L{j+1}"
                match_current, diff_msg_current = compare_lines(
                    ref_lines_orig[i], new_lines_orig[j], line_num_str,
                    file_type, rel_tol, abs_tol
                )

                if match_current:
                    i += 1
                    j += 1
                    continue

                # --- No direct match, try recovery strategies ---
                recovered_this_step = False

                # Strategy 1: Check for a 2-line swap (IGNORE IF MATCHES)
                if i + 1 < len(ref_lines_orig) and j + 1 < len(new_lines_orig):
                    match_r_n1, _ = compare_lines(ref_lines_orig[i], new_lines_orig[j+1], f"Ref L{i+1}/New L{j+2}", file_type, rel_tol, abs_tol)
                    match_r1_n, _ = compare_lines(ref_lines_orig[i+1], new_lines_orig[j], f"Ref L{i+2}/New L{j+1}", file_type, rel_tol, abs_tol)
                    if match_r_n1 and match_r1_n:
                        # This is a 2-line swap, treat as a match for these pairs
                        i += 2
                        j += 2
                        recovered_this_step = True

                if recovered_this_step: continue

                # Strategy 2: Line i in ref is skipped (or extra line j in new) - HARD DIFFERENCE
                # Try to match ref[i+1] with new[j]
                if i + 1 < len(ref_lines_orig):
                    match_r1_n_skip, _ = compare_lines(ref_lines_orig[i+1], new_lines_orig[j], f"Ref L{i+2}/New L{j+1}", file_type, rel_tol, abs_tol)
                    if match_r1_n_skip:
                        differences.append(
                            f"Line Omission/Insertion (Ref L{i+1}): Line present in Ref, not matched in New at current position: '{ref_lines_orig[i].strip()}'"
                        )
                        i += 1
                        recovered_this_step = True

                if recovered_this_step: continue

                # Strategy 3: Line j in new is skipped (or extra line i in ref) - HARD DIFFERENCE
                # Try to match ref[i] with new[j+1]
                if j + 1 < len(new_lines_orig):
                    match_r_n1_skip, _ = compare_lines(ref_lines_orig[i], new_lines_orig[j+1], f"Ref L{i+1}/New L{j+2}", file_type, rel_tol, abs_tol)
                    if match_r_n1_skip:
                        differences.append(
                            f"Line Omission/Insertion (New L{j+1}): Line present in New, not matched in Ref at current position: '{new_lines_orig[j].strip()}'"
                        )
                        j += 1
                        recovered_this_step = True

                if recovered_this_step: continue

                # If no recovery strategy worked, it's a hard mismatch for the original pair
                differences.append(diff_msg_current)
                i += 1
                j += 1

            # Append remaining lines as hard differences
            while i < len(ref_lines_orig):
                if ref_lines_orig[i].strip(): # Only report if non-blank
                    differences.append(f"Extra line at end of Ref file (Ref L{i+1}): '{ref_lines_orig[i].strip()}'")
                i += 1
            while j < len(new_lines_orig):
                if new_lines_orig[j].strip(): # Only report if non-blank
                    differences.append(f"Extra line at end of New file (New L{j+1}): '{new_lines_orig[j].strip()}'")
                j += 1

    except Exception as e:
        print(f"    Error during file comparison for {file_type}: {e}")
        return "error", [str(e)]

    if not differences:
        print(f"    RESULT: {file_type.upper()} files MATCH (2-line swaps ignored).\n")
        return "match", []
    else:
        print(f"    RESULT: {file_type.upper()} files DIFFER.")
        for diff in differences[:10]: # Print up to 10 differences
            print(f"      - {diff}")
        if len(differences) > 10: print(f"      ... and {len(differences) - 10} more differences.\n")
        else: print()
        return "differs", differences


def main_comparison(new_output_dir_name, ref_output_dir_name):
    overall_summary = {"match": 0, "differs": 0,
                       "ref_missing": 0, "new_missing": 0, "error": 0, "total_files": 0}
    examples_with_diffs = []

    suite_base_path = os.getcwd() # or os.path.abspath(os.path.dirname(__file__))
    if not os.path.isdir(suite_base_path):
        print(f"Error: Suite directory '{suite_base_path}' not found.")
        return 2

    print(f"Starting comparison for test suite in: {suite_base_path}")
    print(f"Comparing '{ref_output_dir_name}/' with '{new_output_dir_name}/'")
    print("=" * 70)

    n_comparable_examples = 0

    for category_name in sorted(os.listdir(suite_base_path)):
        #category_path = os.path.join(suite_base_path, category_name)
        category_path = category_name
        if not os.path.isdir(category_path):
            continue

        print(f"\n\n--- Category: {category_name} ---")
        for molecule_basename in sorted(os.listdir(category_path)):
            molecule_example_path = os.path.join(category_path, molecule_basename)
            if not os.path.isdir(molecule_example_path): continue

            print(f"\n  Processing Example: {molecule_basename}")
            example_had_hard_diffs = False

            ref_dir = os.path.join(molecule_example_path, ref_output_dir_name)
            new_dir = os.path.join(molecule_example_path, new_output_dir_name)

            if not os.path.isdir(new_dir):
                print(f"    Skipping {molecule_basename}: New output directory '{new_output_dir_name}' not found.")
                continue
            if not os.path.isdir(ref_dir) : # Reference dir may not exist if no ref files for this example
                print(f"    Note: Reference output directory '{ref_output_dir_name}' not found for {molecule_basename}.")
                continue

            n_comparable_examples += 1
            for ext in [".fit", ".cat", ".out"]:
                overall_summary["total_files"] += 1
                file_type_str = ext[1:] # "fit", "cat", "out"

                ref_file = os.path.join(ref_dir, f"{molecule_basename}{ext}")
                new_file = os.path.join(new_dir, f"{molecule_basename}{ext}")

                # Special handling for .cat and .out which might not always be in reference_outputs
                if ext in [".cat", ".out"] and not os.path.exists(ref_file):
                    print(f"    Reference file {os.path.basename(ref_file)} not found. Skipping comparison for this optional file.")
                    overall_summary["ref_missing"] += 1
                    continue

                status, diff_list = compare_file_pair_final(ref_file, new_file, file_type_str, DEFAULT_REL_TOL, DEFAULT_ABS_TOL)
                overall_summary[status] += 1
                if status == "differs":
                    example_had_hard_diffs = True

            if example_had_hard_diffs:
                examples_with_diffs.append(f"{category_name}/{molecule_basename}")

    print("\n" + "=" * 70)
    print("Overall Comparison Summary:")
    print(f"  Total file pairs processed: {overall_summary['total_files']}")
    print(f"  Matching files:           {overall_summary['match']}")
    print(f"  Differing files:          {overall_summary['differs']}")
    print(f"  Reference files missing:  {overall_summary['ref_missing']}")
    print(f"  New output files missing: {overall_summary['new_missing']}")
    print(f"  Errors during comparison: {overall_summary['error']}")

    if examples_with_diffs:
        print(f"\nExamples with at least one differing file ({len(examples_with_diffs)}/{n_comparable_examples}):")
        for ex in examples_with_diffs:
            print(f"  - {ex}")
    elif overall_summary['differs'] == 0 and overall_summary['new_missing'] == 0 and overall_summary['error'] == 0:
        print("\nAll processed files matched or reference files were appropriately skipped!")
    else:
        print("\nReview summary for details on differences, missing files, or errors.")

    return overall_summary['differs']

if __name__ == "__main__":
    REF_OUTPUT_DIR_NAME = "v2008_results"  # "reference_outputs"

    if len(sys.argv) < 2:
        print("Usage: compare_results.py <comparison_subdir> [baseline_reference_subdir]")
        sys.exit(1)

    sys.exit(main_comparison(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else REF_OUTPUT_DIR_NAME ))