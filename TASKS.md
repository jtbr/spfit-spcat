# SPFIT/SPCAT Modernization Task List

This document outlines the prioritized tasks for modernizing the SPFIT/SPCAT software suite, focusing on enabling programmatic use, improving performance, and enhancing readability/maintainability, while ensuring quantitative results and current functionality are retained.

## Completed Tasks

- [x] **Task 1: Code Readability and Clarity Improvements (Preliminary)**
    - **Description**: Enhance the readability and clarity of the core C code by adding comments, documenting functions, and making targeted, safe variable name changes. This is a preliminary step to facilitate subsequent refactoring.
    - **Implementation Plan**:
        - **Target Files**: `calfit.c`, `dpi.c`, `spinv.c` (split into `spinv_setup.c`, `spinv_spin_symmetry.c`, `spinv_linalg_sort.c`, `spinv_hamiltonian.c`, and `spinv_utils.c`, with common header `spinv_internal.h`).
        - **Actions**:
            - Add inline comments to non-obvious code blocks and complex logic.
            - Document function headers with their purpose, parameters, and return values.
            - Rename variables to be more descriptive, but *only* when their meaning is absolutely clear and the change is low-risk. There can be value in keeping short variable names as well.
        - **Constraint**: Only add comments or change variable names when *certain* of understanding the code in question. Limit changes as much as possible to avoid introducing bugs. DO NOT change the code logic itself.

- [x] **Task 2: Fix Compiler Warnings and Remove Unused Code**
    - **Description**: Address compiler warnings identified during the build process and remove unused functions and files.
    - **Implementation Plan**:
        - **Target Files**: Files with compiler warnings, including `calfit.c`, `subfit.c`, `spinv_setup.c`, `spinv_utils.c`, `ulib.c`, `sortsub.c`, and `calmrg.c`.
        - **Actions**:
            - Remove or properly use variables flagged as "set but not used".
            - Address other warnings such as ignoring return values of functions.
            - Remove unused functions and files identified by static analysis.
        - **Constraint**: Only make minimal changes necessary to fix the warnings without altering the code's behavior.
    - **Results**: All compiler warnings have been fixed, and the code now compiles cleanly with no warnings. Static analysis with `cppcheck` also shows no issues. Unused functions and code files have been removed.


- [x] **Task 3: Obsolete Functions and Portability Issues**
    - **Description**: Addressed remaining code quality issues identified by `cppcheck`, including replacing obsolete functions and fixing portability warnings.
    - **Implementation Plan**:
        - **Obsolete Functions**: All uses of `gets()` have been removed and replaced with safer alternatives like `fgets()`.
        - **Portability Warnings**: Fixed pointer casting issues in `dblas.c` and verified that no other portability issues remain in the codebase.
    - **Results**: The codebase now adheres to modern C standards, with no obsolete functions or unsafe pointer casts. All changes were verified by running the test suite and static analysis.
    - **Constraint**: Quantitative results remain unchanged, and current functionality is retained.

- [x] **Task 6: Build System Modernization**
  - **Description**: Simplify the process of building the software and integrating it with modern development environments.
  - **Implementation Plan**:
    - **Migrate to CMake**: Convert the existing `Makefile` build system to CMake. CMake provides superior cross-platform support, robust dependency tracking, and a more flexible way to define build configurations, including generating project files for various IDEs and creating shared/static libraries.
  - **Constraint**: Quantitative results must remain unchanged, and current functionality must be retained.
  - **Results**: The build system has been migrated to CMake, providing improved cross-platform support and easier integration with modern development environments. Both CMake and Makefile build systems are now supported.

## In Progress Tasks

- [ ] **Task 4: Modularization and Decoupling with C++ Interface**
    - **Description**: Refactor the core SPFIT/SPCAT functionality to be more modular, using C structs/C++ classes for inputs/outputs, and creating a clean C++ interface layer. This is crucial for enabling programmatic use from modern software like Python.
    - **Implementation Plan**: Implementing a Strategy Pattern for Calculation Engines.
        - **Step 1: Create the `CalculationEngine` Interface and Context Structs**:
            - Create the `CalculationEngine.hpp` header file with the abstract base class.
            - Create two new header files, `SpinvContext.hpp` and `DpiContext.hpp`.
            - In `SpinvContext.hpp`, define a `SpinvContext` struct that contains all the global variables used by the `spinv_*.c` files.
            - In `DpiContext.hpp`, define a `DpiContext` struct that contains all the global variables used by `dpi.c`.
        - **Step 2: Create the `SpinvEngine` Class and Refactor `spinv_*.c`**:
            - Create the `SpinvEngine.hpp` and `SpinvEngine.cpp` files.
            - The `SpinvEngine` class will have a member variable of type `SpinvContext`.
            - Modify the functions in the `spinv_*.c` files to take a pointer to the `SpinvContext` struct as a parameter, instead of accessing global variables directly.
            - The `SpinvEngine`'s methods will then call these modified C functions, passing in the `SpinvContext` member variable.
        - **Step 3: Create the `DpiEngine` Class and Refactor `dpi.c`**:
            - Create the `DpiEngine.hpp` and `DpiEngine.cpp` files.
            - The `DpiEngine` class will have a member variable of type `DpiContext`.
            - Modify the functions in `dpi.c` to take a pointer to the `DpiContext` struct as a parameter.
            - The `DpiEngine`'s methods will then call these modified C functions, passing in the `DpiContext` member variable.
        - **Step 4: Integrate the `CalculationEngine` into `spfit` and `spcat`**:
            - Modify the `main` functions in `calfit.c` and `calcat.c` to use the `CalculationEngine` interface.
        - **Step 5: Preliminary simplifications**:
            - Remove `zero` and `szero` from `SpinvContext` and `DpiContext`, use `memset` instead of `dcopy()` where simple, remove context parameter from functions that don't need it.
        - **Step 6: Refactor `fit_main.cpp` and `cat_main.cpp`**:
            - Refactor `fit_main.cpp` and `cat_main.cpp` to be more modular.
            - This will involve creating `Spfit` and `Spcat` classes and moving the logic from the `main` functions into these classes.
            In this step we'll need to:
              - **Define Clear C++ APIs**: Create well-defined C++ class interfaces (e.g., `Spfit`, `Spcat`) with methods that encapsulate core functionalities. These methods should take C++ data structures as input and return C++ data structures as output.
              - **Separate I/O from Logic**: Modify core calculation functions to operate on in-memory data structures (passed via the new C++/C structs/classes) rather than directly performing file I/O. Create separate C++ utility functions or methods for reading from and writing to files.
              - **Leverage C++ Features**: Utilize RAII (Resource Acquisition Is Initialization) with `std::unique_ptr` and `std::vector` for robust memory management within the new C++ interface layer.
    - **Progress**: Steps 1-5 are complete. Step 6 is partially complete with significant progress on `CalFit` class (renamed from `Spfit` for clarity).
      - **Completed in Step 6**:
        - Created `CalFit` class with methods for parameter initialization (`initializeParameters`), spectral line processing (`processLines`), and formatting (`parer`, `qnfmt2`).
        - Implemented core logic for reading and processing spectral lines (`linein`, `lineix`, `getblk`).
        - Updated `fit_main.cpp` to use the `CalFit` class integrated with the calculation engine.
        - Implemented iterative fitting logic in `performIteration` using the Marquardt-Levenberg algorithm.
        - Implemented error calculation and statistics in `calculateErrors`.
        - Completed I/O separation by implementing `readInput` and `writeOutput` in `CalFitIO` for handling file operations (these are not fully utilized yet; there are still scratch files being used.)
      - **Remaining in Step 6**:
        - Debug current issues in reproducing correct functionality (a number of discrepancies in output files remain, compared with the baseline outputs.)
        - Once this is validated, commit a new baseline to git
        - Refactor `cat_main.cpp` into a `CalCat` class following the same modularization pattern.
        - Make full use of `CalFitIO` and remove use of scratch files.
        - Note that some input/output files are shared between `fit` and `cal` and it may make sense that IO is modularized by file type and shared between `CalCat` and `CalFit`. This needs to be explored first.
        - we may need to refactor linein and getlin and possibly other code to remove dependencies on file IO and separate that from the logic
    - **Details**: The refactoring has focused on disentangling file I/O from calculation logic within the `CalFit` class. Temporary file (or in-memory file) placeholders are in use until full I/O separation is achieved.
    - **Plan**:
        1.  **Spfit/Spcat Class Structure:**
            *   Create `Spfit` (implemented as `CalFit`) and `Spcat` classes.
            *   These classes will encapsulate the core fitting and catalog generation logic, respectively.
            *   Each class will have a `run` method that performs the main calculation.

        2.  **Input/Output Data Structures:**
            *   Define C++ structs or classes to hold input parameters and output results. These structures should be compatible with potential Python bindings (e.g., using standard data types).
            *   These data structures will be used to pass data between the `Spfit`/`Spcat` classes and the calling code.
            *   This will allow the software to be run without files, by passing data directly to the `run` method.

        3.  **Separation of I/O and Calculation Logic:**
            *   Create separate functions or classes for reading input data from files and writing output data to files.
            *   The `Spfit`/`Spcat` classes will operate on in-memory data structures, not directly on files.
            *   This will involve modifying the `SpinvEngine`, `DpiEngine`, `spinv_*.c`, and `dpi.c` files to separate I/O from calculation.

        4.  **C++ Interface Design:**
            *   Design a clear and simple C++ API for the `Spfit` and `Spcat` classes.
            *   The API should expose all current functionality, but with sensible defaults to make it easy to use.
            *   The API should allow the user to specify input parameters and output options.

        5.  **File-Based Execution:**
            *   Maintain the existing ability to run `spfit` and `spcat` using a set of input files.
            *   The `main` function will parse command-line arguments and call the appropriate I/O functions to read input data from files.
            *   The input data will then be passed to the `Spfit`/`Spcat` classes for calculation.

        6.  **Amortization of Startup Costs:**
            *   Design the `Spfit` and `Spcat` classes to minimize startup costs.
            *   This may involve caching data structures or pre-calculating values that are used repeatedly.
            *   This will allow the software to be used efficiently in situations where it is called multiple times.

        7.  **Testing and Validation:**
            *   Leverage the existing test suite (input files and expected outputs) for validation.
            *   Add unit tests as needed to cover specific functionality.
    - **Expected Outcome**: A highly modular codebase with a clean C++ API, enabling easy programmatic integration and future Python bindings.
    - **Constraint**: Quantitative results must remain unchanged, and current functionality must be retained.

- [ ] **Task 5: Performance Optimizations**
    - **Description**: Speed up calculations and exploit modern hardware.
    - **Implementation Plan**:
        - **Profile Critical Sections**: Use profiling tools (e.g., `gprof` or `perf`) to identify the most time-consuming parts of the code, especially matrix operations (Hamiltonian diagonalization, least-squares fitting).
        - **Parallelization**: Explore opportunities for parallelizing independent calculations, particularly within matrix operations. OpenMP could be a suitable choice for C.
        - **Leverage Optimized Libraries**: Ensure optimal and correct use of external BLAS/LAPACK libraries (like OpenBLAS) for all linear algebra operations, as these are highly optimized.
        - **Memory Access Patterns**: Review and optimize data structures and access patterns to improve cache utilization and reduce memory access latency.
    - **Constraint**: Quantitative results must remain unchanged, and current functionality must be retained.


- [ ] **Task 7: Error Handling and Input Validation**
    - **Description**: Improve robustness and provide better diagnostics for programmatic users.
    - **Implementation Plan**:
        - **Return Error Codes**: Modify functions to return explicit error codes or status indicators instead of calling `exit()`, allowing the calling program/script to handle errors gracefully.
        - **Informative Logging**: Implement a simple, configurable logging mechanism to output detailed error messages, warnings, and debug information to `stderr` or a log file, providing context for issues.
        - **Input Validation**: Add comprehensive checks at function boundaries to validate input parameters and data structures, preventing crashes due to malformed data.
        - Consider adding more user-friendly features or command line arguments
    - **Constraint**: Quantitative results must remain unchanged, and current functionality must be retained.

- [ ] **Task 8: Memory Management Refinement**
    - **Description**: Improve memory management practices to reduce leaks and improve stability.
    - **Implementation Plan**:
        - As part of the modularization effort (Task 3), ensure clear ownership and lifecycle management of dynamically allocated memory. This might involve designing APIs where callers provide buffers, or where modules explicitly manage their internal memory and provide cleanup functions.
    - **Constraint**: Quantitative results must remain unchanged, and current functionality must be retained.

- [ ] - **Task 9: Python interface and example usage**
    - Create an example for using the C++ interface from C++, showcasing important features and optimizations for using it
    - Create a python interface based upon the C++ interface (`pybind11`),
    - Create example code for setting up and using the library from python

- [ ] **Task 10: Documentation Improvement**:
    - Enhance code documentation
    - Create better user documentation
    - Document the scientific principles and algorithms

## Implementation Plan

The tasks will be executed sequentially, with thorough testing and verification after each major change to ensure the core requirements are met.

### Relevant Files

*   `TASKS.md` - This task list itself.
*   `memory-bank/activeContext.md` - Will be updated with progress and insights.
*   `memory-bank/progress.md` - Will be updated with overall project status.
*   Core C files: `calfit.c`, `dpi.c`, `spinv_setup.c`, `spinv_spin_symmetry.c`, `spinv_linalg_sort.c`, `spinv_hamiltonian.c`, `spinv_utils.c`, `spinv_internal.h` (related to Task 1).
*   Files identified by `cppcheck` for issues or potential removal (for Task 2).
*   New C++ header and source files: `CalculationEngine.hpp`, `SpinvEngine.hpp`, `SpinvEngine.cpp`, `SpinvContext.hpp`, `DpiEngine.hpp`, `DpiEngine.cpp`, `DpiContext.hpp` (related to Task 4).
