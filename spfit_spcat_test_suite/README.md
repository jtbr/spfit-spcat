# SPFIT/SPCAT Test Suite

This directory contains a collection of test cases for the SPFIT and SPCAT programs,
based on examples from the [CDMS Pickett software examples page](https://cdms.astro.uni-koeln.de/classic/predictions/pickett/beispiele/).

## Purpose

The purpose of this suite is to:
1. Provide a standardized set of inputs for testing SPFIT and SPCAT installations.
2. Allow verification of results by comparing generated outputs against reference outputs.

## Directory Structure

The test suite is organized as follows:

spfit_spcat_test_suite/
├── README.md
├── [category_name]/ # e.g., diatomic_molecules
│ ├── [molecule_example_name]/ # e.g., CO
│ │ ├── [basename].par # Input parameter file for SPFIT
│ │ ├── [basename].lin # Input line file for SPFIT
│ │ ├── [basename].int # Input intensity file for SPCAT
│ │ └── reference_outputs/
│ │ ├── [basename].fit # Reference output from SPFIT
│ │ ├── [basename].cat # Reference output from SPCAT, if available
│ │ └── [basename].out # Reference general output from SPCAT, if available
│ │ └── v2008_outputs/
│ │ ├── [basename].fit # Reference output from SPFIT
│ │ ├── [basename].cat # Reference output from SPCAT
│ │ └── [basename].out # Reference general output from SPCAT

Note: I have confirmed that differences between the reference and v2008 versions don't appear to be due to compiler flags or due to changes from laser_kelvin.

## Running Tests

**Prerequisites:**
* SPFIT and SPCAT programs pre-built in code directory

Navigate to the `spfit_spcat_test_suite` subdirectory and run the `run_tests` script.
```sh
cd spfit_spcat_test_suite && python run_tests.py [output_subdir]
```

This will run spfit and spcat on copies of the input files for each example in each example category, and put the results in the `output_subdir` subdirectory of each example directory. It runs from a temporary directory to avoid overwriting the input files.

## Verifying Results

After running SPFIT and SPCAT, compare the generated output files (e.g., `[basename].fit`, `[basename].cat`) with the files in a reference subdirectory (e.g. `reference_outputs/`)

*   **Text-based comparison (for `.fit`, `.out`):**
    ```bash
    diff [basename].fit reference_outputs/[basename].fit
    diff [basename].out reference_outputs/[basename].out
    diff [basename].cat reference_outputs/[basename].cat
    ```
    (I recommend a GUI diff tool for these files, as there may be very small changes)

* **Automated analysis of all example outputs**
    ```bash
    python compare_results.py [reference_subdir] [comparison_subdir]
    ```

    This will produce a report for each example and each output file, detailing differences and noting matches.

## Notes

*   SPFIT and SPCAT overwrite existing output files in the current directory. Each test case is in its own directory to prevent interference.
*   The reference output files were downloaded from the CDMS Pickett examples website.
