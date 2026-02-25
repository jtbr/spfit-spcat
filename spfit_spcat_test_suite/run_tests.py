import os
import subprocess
import shutil
import tempfile
import time
import sys

# spcat/spfit executables are in the parent directory to this directory
EXE_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
spfit_path = os.path.join(EXE_PATH, "spfit")
spcat_path = os.path.join(EXE_PATH, "spcat")
#OUTPUT_SUBDIR_NAME = "v2008_outputs" # Name of the subdirectory to store outputs

# Define expected output file extensions
# SPFIT also overwrites the .par file and creates a .bak from the original .par
SPFIT_OUTPUT_EXTENSIONS = [".fit", ".var", ".bak", ".par"]
SPCAT_OUTPUT_EXTENSIONS = [".cat", ".out", ".egy", ".str"] # .egy and .str are optional

def run_commands_and_move_outputs_dynamic(output_subdir_name, suite_base_path='.', skip=[]):
    original_cwd = os.getcwd()

    # abort if the executables aren't there
    if not os.path.exists(spfit_path):
      print("Error: unable to find spfit executable: ", spfit_path)
      return

    print(f"Starting processing (Dynamic Discovery) for test suite in {suite_base_path}")
    print(f"Outputs will be placed in '{output_subdir_name}' subdirectories for each example.")
    print("-" * 50)

    # Dynamically discover categories (subdirectories of suite_base_path)
    for category_name in sorted(os.listdir(suite_base_path)):
        category_path = os.path.join(suite_base_path, category_name)
        if not os.path.isdir(category_path):
            continue # Skip files like README.md in the root

        print(f"\nProcessing Category: {category_name}")

        # Dynamically discover molecule examples (subdirectories of category_path)
        for local_basename in sorted(os.listdir(category_path)):
            molecule_example_path = os.path.join(category_path, local_basename)
            molecule_example_abspath = os.path.abspath(molecule_example_path)

            if not os.path.isdir(molecule_example_path):
                continue # Skip any files within category directories

            if local_basename in skip:
                print(f"  SKIPPING Molecule: {local_basename} in skip list (--all to disable)")
                continue

            print(f"  Processing Molecule: {local_basename} in {molecule_example_path}")

            # --- Change to the molecule's directory ---
            try:
                # copy input files to a temp dir, as they will be overwritten
                with tempfile.TemporaryDirectory() as temp_dir:
                    # Copy all files from molecule_example_path to temp_dir
                    for filename in os.listdir(molecule_example_path):
                        src_path = os.path.join(molecule_example_path, filename)
                        if os.path.isfile(src_path):
                            dst_path = os.path.join(temp_dir, filename)
                            shutil.copy2(src_path, dst_path)
                            os.chmod(dst_path, 0o644)

                    os.chdir(temp_dir)
                    #print(f"    Changed CWD to: {os.getcwd()}")

                    # output_dir_path is inside the current molecule_example_path
                    output_dir_path_abs = os.path.join(molecule_example_abspath, output_subdir_name)
                    os.makedirs(output_dir_path_abs, exist_ok=True)

                    # --- Run SPFIT ---
                    spfit_success = False
                    spfit_command = [spfit_path, local_basename]
                    print(f"    Running SPFIT: {' '.join(spfit_command)}")
                    try:
                        # Ensure input .par file exists for spfit
                        if not os.path.exists(f"{local_basename}.par"):
                            print(f"    ERROR: SPFIT input file {local_basename}.par not found in {os.getcwd()}. Skipping SPFIT/SPCAT for this molecule.")
                            os.chdir(original_cwd) # Change back before continuing to next molecule
                            continue

                        spfit_log_path = os.path.join(output_dir_path_abs, f"{local_basename}_spfit.log")
                        spfit_result = subprocess.run(spfit_command, capture_output=True, text=True, check=False, timeout=600)
                        with open(spfit_log_path, "w") as f:
                            f.write("--- SPFIT STDOUT ---\n")
                            f.write(spfit_result.stdout)
                            f.write("\n--- SPFIT STDERR ---\n")
                            f.write(spfit_result.stderr)
                        if spfit_result.returncode == 0:
                            print("    SPFIT completed successfully.")
                            spfit_success = True
                        else:
                            print(f"    SPFIT failed with return code {spfit_result.returncode}.")
                            print(f"    SPFIT output logged to {spfit_log_path}")
                    except subprocess.TimeoutExpired:
                        print("    SPFIT command timed out.")
                    except Exception as e:
                        print(f"    An error occurred running SPFIT: {e}")

                    # --- Run SPCAT (only if SPFIT was successful) ---
                    spcat_success = False
                    if spfit_success:
                        if not os.path.exists(f"{local_basename}.var"):
                            print(f"    ERROR: SPCAT input file {local_basename}.var not found (expected from SPFIT). Skipping SPCAT.")
                        elif not os.path.exists(f"{local_basename}.int"):
                            print(f"    ERROR: SPCAT input file {local_basename}.int not found. Skipping SPCAT.")
                        else:
                            spcat_command = [spcat_path, local_basename]
                            print(f"    Running SPCAT: {' '.join(spcat_command)}")
                            try:
                                spcat_log_path = os.path.join(output_dir_path_abs, f"{local_basename}_spcat.log")
                                spcat_result = subprocess.run(spcat_command, capture_output=True, text=True, check=False, timeout=1000)
                                with open(spcat_log_path, "w") as f:
                                    f.write("--- SPCAT STDOUT ---\n")
                                    f.write(spcat_result.stdout)
                                    f.write("\n--- SPCAT STDERR ---\n")
                                    f.write(spcat_result.stderr)
                                if spcat_result.returncode == 0:
                                    print("    SPCAT completed successfully.")
                                    spcat_success = True
                                else:
                                    print(f"    SPCAT failed with return code {spcat_result.returncode}.")
                                    print(f"    SPCAT output logged to {spcat_log_path}")
                            except subprocess.TimeoutExpired:
                                print("    SPCAT command timed out.")
                            except Exception as e:
                                print(f"    An error occurred running SPCAT: {e}")
                    else:
                        print("    Skipping SPCAT due to SPFIT failure or non-completion.")

                    # --- Move all relevant output files ---
                    files_to_move = []
                    if spfit_success: # If SPFIT at least ran, try to move its outputs
                        for ext in SPFIT_OUTPUT_EXTENSIONS:
                            files_to_move.append(f"{local_basename}{ext}")

                    if spcat_success: # Only try to move spcat outputs if it ran successfully
                        for ext in SPCAT_OUTPUT_EXTENSIONS:
                            files_to_move.append(f"{local_basename}{ext}")

                    files_to_move = sorted(list(set(files_to_move))) # Ensure uniqueness

                    for output_file_basename in files_to_move:
                        # output_file_basename is e.g. "co.fit"
                        # source is in current CWD (molecule_example_path)
                        source_file_path = os.path.join(temp_dir, output_file_basename)

                        if os.path.exists(source_file_path):
                            try:
                                destination_file_path = os.path.join(output_dir_path_abs, output_file_basename)
                                shutil.move(source_file_path, destination_file_path)
                                print(f"      Moved {output_file_basename} to {destination_file_path}")
                            except Exception as e:
                                print(f"      Error moving {output_file_basename} to {output_dir_path_abs}: {e}")
                        else:
                            is_optional_spcat_output = spcat_success and \
                                                      any(output_file_basename.endswith(ext) for ext in [".egy", ".str"])
                            if not is_optional_spcat_output: # Don't warn for missing optional files if spcat ran
                                print(f"      Output file {output_file_basename} not found for moving (may be normal if command failed or file is optional).")

            except Exception as e_outer:
                print(f"  Outer error processing {local_basename}: {e_outer}")
            finally:
                os.chdir(original_cwd) # Crucially change back to original CWD
            time.sleep(0.1)

        print("-" * 50)
    print("\nAll processing finished.")

# TODO: allow running one or more specific molecules, to be specified after subdirectory

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: run_tests.py [--all] <subdirectory>")
        print(f"Where the required subdirectory is the one to use for each example's outputs (placed into each molecule directory)")
        print(f"      and --all enables the `clclo2` molecule test, which is very slow")
        sys.exit(1)

    # NOTE: Skip clclo2 example by default as it's too slow to iterate; use only for full tests
    subdir = sys.argv[1]
    skip = ['clclo2']
    if subdir == "--all":
        skip = []
        subdir = sys.argv[2]
    run_commands_and_move_outputs_dynamic(subdir, skip=skip)