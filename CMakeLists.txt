cmake_minimum_required(VERSION 3.5)
project(spfit LANGUAGES C CXX DESCRIPTION "SPFIT/SPCAT")
cmake_policy(SET CMP0057 NEW)
# NOTE: Unlike the hand-maintained Makefile, CMake's `make` and `make all` are identical (both build the
# default ALL target). There is no cmake equivalent to the Makefile's separate `make` vs `make all` targets.
# To build legacy utilities via cmake, either use `cmake -DBUILD_ALL=ON` or `make <target>` individually.
# Use the root Makefile directly if you need `make all` semantics (builds all executables by default).
option(BUILD_ALL "Build all executables including legacy utilities (cmake -DBUILD_ALL=ON)" OFF)

# TODO re-add debug options
# NOTE: -ffp-contract=off disables fused multiply-add (FMA), which change the accumulation order of BLAS dot
# products, leading to sign-flipped columns in the packed Cholesky factor. This is not essential (results are
# physically equivalent), but together with use of dblas.c allows exact matching with v2008_baseline.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -march=native -ffp-contract=off")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O3 -march=native -ffp-contract=off")

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

# Find BLAS
option(USE_SYSTEM_BLAS "Use system BLAS library (OpenBLAS, etc.)" OFF)
if(USE_SYSTEM_BLAS)
  find_package(BLAS) # Removed REQUIRED to allow fallback
else()
  set(BLAS_FOUND FALSE)
endif()
if (BLAS_FOUND)
    message(STATUS "BLAS found, using system BLAS: ${BLAS_LIBRARIES}")
    include_directories(${BLAS_INCLUDE_DIRS}) # Only include if BLAS is found
else()
    message(STATUS "BLAS not found. Relevant targets will compile dblas.c as a fallback.")
endif()

# Define all potential source files
# This list is used for set_source_files_properties
set(SOURCE_FILES
    src/spfit/fit_main.cpp
    src/spfit/subfit.cpp
    src/spfit/CalFit.cpp
    src/spfit/CalFitIO.cpp
    src/spfit/CalFit_helpers.cpp
    src/spcat/cat_main.cpp
    src/spcat/sortsub.c
    src/spcat/CalCat.cpp
    src/spcat/CalCatIO.cpp
    src/spcat/CalCat_helpers.cpp
    src/engine/spinv_setup.cpp
    src/engine/spinv_spin_symmetry.cpp
    src/engine/spinv_linalg_sort.c
    src/engine/spinv_hamiltonian.cpp
    src/engine/spinv_utils.cpp
    src/engine/spinit.cpp
    src/engine/SpinvEngine.cpp
    src/engine/DpiEngine.cpp
    src/engine/dpi.cpp
    src/splib/ulib.c
    src/splib/cnjj.c
    src/splib/lsqfit.c
    src/splib/dblas.c
    src/splib/ftran.c
    src/splib/catutil.c
    src/common/Logger.cpp
    src/common/SigintFlag.cpp
    src/common/file_helpers.cpp
    src/legacy_apps/calmrg.cpp
    src/legacy_apps/calbak.cpp
    src/legacy_apps/termval.cpp
    src/legacy_apps/stark.cpp
    src/legacy_apps/moiam.cpp
    src/legacy_apps/iamcalc.cpp
    src/legacy_apps/sortn.cpp
    src/legacy_apps/reassign.cpp
    src/legacy_apps/sortegy.cpp
    src/legacy_apps/iambak.cpp
    src/legacy_apps/readopt.c
)

set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers")
#set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "-Wno-unused-function")

# Add executable macro
macro(add_executable_plus name uses_blas_fallback)
  # The actual source files for the target are now in ARGN
  set(current_sources ${ARGN})

  if(NOT BLAS_FOUND AND ${uses_blas_fallback})
    list(APPEND current_sources src/splib/dblas.c) # Add dblas.c if no system BLAS and target needs it
    message(STATUS "Target ${name}: System BLAS not found. Compiling dblas.c as fallback.")
  elseif(BLAS_FOUND AND ${uses_blas_fallback})
    message(STATUS "Target ${name}: Using system BLAS.")
  elseif(NOT ${uses_blas_fallback})
    message(STATUS "Target ${name}: Does not require BLAS.")
  endif()

  add_executable(${name} ${current_sources})

  # Link system BLAS if found and the target is meant to use BLAS
  if(BLAS_FOUND AND ${uses_blas_fallback} AND BLAS_LIBRARIES)
    target_link_libraries(${name} PRIVATE ${BLAS_LIBRARIES})
  endif()

  if ("${name}" MATCHES "^(spfit|spcat)$")
    target_link_libraries(${name} PRIVATE m stdc++)
  else()
    target_link_libraries(${name} PRIVATE m)
  endif()
  set_property(TARGET ${name} PROPERTY LINKER_LANGUAGE CXX)
endmacro()

# Default build: spfit spcat calmrg. Legacy utilities excluded unless -DBUILD_ALL=ON.
set(OPTIONAL_EXECUTABLES moiam stark termval sortn calbak reassign sortegy iambak iamcalc)

set(ENGINE_SRCS
    src/engine/spinv_setup.cpp
    src/engine/spinv_spin_symmetry.cpp
    src/engine/spinv_linalg_sort.c
    src/engine/spinv_hamiltonian.cpp
    src/engine/spinv_utils.cpp
    src/engine/spinit.cpp
    src/engine/dpi.cpp
    src/engine/SpinvEngine.cpp
    src/engine/DpiEngine.cpp
)
set(SPLIB_SRCS
    src/splib/catutil.c
    src/splib/ulib.c
    src/splib/cnjj.c
    src/splib/lsqfit.c
)
set(COMMON_SRCS
    src/common/Logger.cpp
    src/common/SigintFlag.cpp
    src/common/file_helpers.cpp
)

# Add spfit
add_executable_plus(spfit TRUE
    src/spfit/fit_main.cpp
    src/spfit/CalFit.cpp
    src/spfit/CalFit_helpers.cpp
    src/spfit/CalFitIO.cpp
    src/spfit/subfit.cpp
    src/splib/ftran.c
    ${ENGINE_SRCS}
    ${SPLIB_SRCS}
    ${COMMON_SRCS}
)

# Add spcat
add_executable_plus(spcat TRUE
    src/spcat/cat_main.cpp
    src/spcat/CalCat.cpp
    src/spcat/CalCat_helpers.cpp
    src/spcat/CalCatIO.cpp
    src/spcat/sortsub.c
    ${ENGINE_SRCS}
    ${SPLIB_SRCS}
    ${COMMON_SRCS}
)

# Shared helpers required by all legacy executables
set(CPP_HELPERS ${COMMON_SRCS})
set(SPLIB ${SPLIB_SRCS})

# Default: calmrg (matches Makefile default targets: spfit spcat calmrg)
add_executable_plus(calmrg TRUE
    src/legacy_apps/calmrg.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)

# Optional executables — always defined so `make <target>` works individually;
# excluded from the default `make all` unless -DBUILD_ALL=ON.
add_executable_plus(calbak TRUE
    src/legacy_apps/calbak.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(termval TRUE
    src/legacy_apps/termval.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(stark TRUE
    src/legacy_apps/stark.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(moiam TRUE
    src/legacy_apps/moiam.cpp
    src/splib/ftran.c
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(iamcalc TRUE
    src/legacy_apps/iamcalc.cpp
    src/splib/ftran.c
    src/legacy_apps/readopt.c
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(sortn FALSE
    src/legacy_apps/sortn.cpp
    src/spcat/sortsub.c
    ${CPP_HELPERS}
)
add_executable_plus(reassign TRUE
    src/legacy_apps/reassign.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(sortegy TRUE
    src/legacy_apps/sortegy.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(iambak TRUE
    src/legacy_apps/iambak.cpp
    src/legacy_apps/readopt.c
    ${CPP_HELPERS}
    ${SPLIB}
)

# Exclude legacy utilities from default build unless -DBUILD_ALL=ON.
# Note: targets are still defined above, so `make <target>` always works individually.
if(NOT BUILD_ALL)
  foreach(target ${OPTIONAL_EXECUTABLES})
    set_target_properties(${target} PROPERTIES EXCLUDE_FROM_ALL TRUE)
  endforeach()
endif()
