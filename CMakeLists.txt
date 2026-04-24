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

include_directories(.)

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
    fit_main.cpp
    subfit.cpp
    spinv_setup.cpp
    spinv_spin_symmetry.cpp
    spinv_linalg_sort.c
    spinv_hamiltonian.cpp
    spinv_utils.cpp
    spinit.cpp
    CalFit.cpp
    CalFitIO.cpp
    SpinvEngine.cpp
    DpiEngine.cpp
    cat_main.cpp
    sortsub.c
    calmrg.cpp
    calbak.cpp
    termval.cpp
    stark.cpp
    moiam.cpp
    iamcalc.cpp
    sortn.cpp
    reassign.cpp
    sortegy.cpp
    iambak.cpp
    ulib.c
    cnjj.c
    lsqfit.c
    dblas.c    # dblas.c is in this list
    dpi.cpp
    ftran.c
    catutil.c
    readopt.c
    Logger.cpp
    SigintFlag.cpp
    file_helpers.cpp
)

set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "-Wno-missing-field-initializers")
#set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "-Wno-unused-function")

# Add executable macro
macro(add_executable_plus name uses_blas_fallback)
  # The actual source files for the target are now in ARGN
  set(current_sources ${ARGN})

  if(NOT BLAS_FOUND AND ${uses_blas_fallback})
    list(APPEND current_sources dblas.c) # Add dblas.c if no system BLAS and target needs it
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

# Add spfit
if(TRUE)
  add_executable_plus(spfit TRUE # This target uses BLAS
    fit_main.cpp
    CalFit.cpp
    CalFit_helpers.cpp
    CalFitIO.cpp
    subfit.cpp
    spinv_setup.cpp
    spinv_spin_symmetry.cpp
    spinv_linalg_sort.c
    spinv_hamiltonian.cpp
    spinv_utils.cpp
    spinit.cpp
    dpi.cpp
    SpinvEngine.cpp
    DpiEngine.cpp
    catutil.c
    ulib.c
    cnjj.c
    lsqfit.c
    # dblas.c removed - macro handles it
    ftran.c
    Logger.cpp
    SigintFlag.cpp
    file_helpers.cpp
  )
endif()

# Add spcat
if(TRUE)
  add_executable_plus(spcat TRUE # This target uses BLAS
    cat_main.cpp
    CalCat.cpp
    CalCat_helpers.cpp
    CalCatIO.cpp
    sortsub.c
    spinv_setup.cpp
    spinv_spin_symmetry.cpp
    spinv_linalg_sort.c
    spinv_hamiltonian.cpp
    spinv_utils.cpp
    spinit.cpp
    dpi.cpp
    SpinvEngine.cpp
    DpiEngine.cpp
    catutil.c
    ulib.c
    cnjj.c
    lsqfit.c
    # dblas.c removed - macro handles it
    Logger.cpp
    SigintFlag.cpp
    file_helpers.cpp
  )
endif()

# # Add dpfit
# if (NOT DEFINED EXCLUDE_CALMRG AND "calmrg" IN_LIST EXECUTABLES)
#   add_executable_plus(calmrg TRUE # This target uses BLAS
#       calfit.c
#       sortfit.c
#       dpi.c
#       ulib.c cnjj.c slibgcc.c catutil.c lsqfit.c # splib
#       # dblas.c removed
#   )
# endif()

# # Add dpcat
# if (NOT DEFINED EXCLUDE_CALMRG AND "calmrg" IN_LIST EXECUTABLES)
#   add_executable_plus(calmrg TRUE # This target uses BLAS
#       calcat.c
#       sortsub.c
#       dpi.c
#       ulib.c cnjj.c slibgcc.c catutil.c lsqfit.c # splib
#       # dblas.c removed
#   )
# endif()

# Shared helpers required by all legacy executables
set(CPP_HELPERS file_helpers.cpp SigintFlag.cpp Logger.cpp)
set(SPLIB ulib.c cnjj.c catutil.c lsqfit.c)

# Default: calmrg (matches Makefile default targets: spfit spcat calmrg)
add_executable_plus(calmrg TRUE
    calmrg.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)

# Optional executables — always defined so `make <target>` works individually;
# excluded from the default `make all` unless -DBUILD_ALL=ON.
add_executable_plus(calbak TRUE
    calbak.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(termval TRUE
    termval.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(stark TRUE
    stark.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(moiam TRUE
    moiam.cpp
    ftran.c
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(iamcalc TRUE
    iamcalc.cpp
    ftran.c
    readopt.c
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(sortn FALSE
    sortn.cpp
    sortsub.c
    ${CPP_HELPERS}
)
add_executable_plus(reassign TRUE
    reassign.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(sortegy TRUE
    sortegy.cpp
    ${CPP_HELPERS}
    ${SPLIB}
)
add_executable_plus(iambak TRUE
    iambak.cpp
    readopt.c
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
