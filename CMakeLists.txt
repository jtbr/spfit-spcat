cmake_minimum_required(VERSION 3.0)
project(spfit LANGUAGES C CXX DESCRIPTION "SPFIT/SPCAT")
cmake_policy(SET CMP0057 NEW)
# set with `cmake [path] -DBUILD_ALL=ON`
option(BUILD_ALL "Build all executables" OFF)

# TODO re-add debug options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -march=native")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O3 -march=native")

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

# Find BLAS
find_package(BLAS) # Removed REQUIRED to allow fallback
if (BLAS_FOUND)
    message(STATUS "BLAS found, using system BLAS: ${BLAS_LIBRARIES}")
    include_directories(${BLAS_INCLUDE_DIRS}) # Only include if BLAS is found
else()
    message(STATUS "BLAS not found. Relevant targets will compile dblas.c as a fallback.")
endif()

# Define all potential source files
# This list is used for set_source_files_properties
set(SOURCE_FILES
    calfit.cpp
    subfit.c
    spinv_setup.c
    spinv_spin_symmetry.c
    spinv_linalg_sort.c
    spinv_hamiltonian.c
    spinv_utils.c
    spinit.c
    SpinvEngine.cpp
    calcat.cpp
    sortsub.c
    calmrg.c
    calbak.c
    termval.c
    stark.c
    moiam.c
    iamcalc.c  # Even if target is commented out, keep for properties if file exists
    cnvwn.c
    sortn.c
    reassign.c
    sortegy.c
    iambak.c
    ulib.c
    cnjj.c
    slibgcc.c
    lsqfit.c
    dblas.c    # dblas.c is in this list
    dpi.c
    ftran.c
    catutil.c  # Added missing common files
    readopt.c  # Added missing common files
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

# Define executables
set(DEFAULT_EXECUTABLES spfit spcat)
set(OPTIONAL_EXECUTABLES calmrg moiam stark termval sortn calbak reassign sortegy iambak)
if(BUILD_ALL)
    set(EXECUTABLES ${DEFAULT_EXECUTABLES} ${OPTIONAL_EXECUTABLES})
else()
    set(EXECUTABLES ${DEFAULT_EXECUTABLES})
endif()


if(NOT DEFINED EXECUTABLES)
  set(EXECUTABLES ${DEFAULT_EXECUTABLES})
endif()

# Add spfit
if("spfit" IN_LIST EXECUTABLES)
  add_executable_plus(spfit TRUE # This target uses BLAS
    calfit.cpp
    subfit.c
    spinv_setup.c
    spinv_spin_symmetry.c
    spinv_linalg_sort.c
    spinv_hamiltonian.c
    spinv_utils.c
    spinit.c
    dpi.c
    SpinvEngine.cpp
    DpiEngine.cpp
    catutil.c
    ulib.c
    cnjj.c
    slibgcc.c
    lsqfit.c
    # dblas.c removed - macro handles it
    ftran.c
  )
endif()

# Add spcat
if("spcat" IN_LIST EXECUTABLES)
  add_executable_plus(spcat TRUE # This target uses BLAS
    calcat.cpp
    sortsub.c
    spinv_setup.c
    spinv_spin_symmetry.c
    spinv_linalg_sort.c
    spinv_hamiltonian.c
    spinv_utils.c
    spinit.c
    dpi.c
    SpinvEngine.cpp
    DpiEngine.cpp
    catutil.c
    ulib.c
    cnjj.c
    slibgcc.c
    lsqfit.c
    # dblas.c removed - macro handles it
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

# --non defaults --

# Add calmrg
if (NOT DEFINED EXCLUDE_CALMRG AND "calmrg" IN_LIST EXECUTABLES)
  add_executable_plus(calmrg TRUE # This target uses BLAS
      calmrg.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      # dblas.c removed
  )
endif()

# Add calbak
if (NOT DEFINED EXCLUDE_CALBAK AND "calbak" IN_LIST EXECUTABLES)
  add_executable_plus(calbak TRUE # This target uses BLAS
      calbak.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      # dblas.c removed
  )
endif()

# Add termval
if (NOT DEFINED EXCLUDE_TERMVAL AND "termval" IN_LIST EXECUTABLES)
  add_executable_plus(termval TRUE # This target uses BLAS
      termval.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      # dblas.c removed
  )
endif()

# Add stark
if (NOT DEFINED EXCLUDE_STARK AND "stark" IN_LIST EXECUTABLES)
  add_executable_plus(stark TRUE # This target uses BLAS
      stark.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      # dblas.c removed
  )
endif()

# Add moiam
if (NOT DEFINED EXCLUDE_MOIAM AND "moiam" IN_LIST EXECUTABLES)
  add_executable_plus(moiam TRUE # This target uses BLAS
      moiam.c
      ftran.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      # dblas.c removed
  )
endif()

# # Add iamcalc (NOT WORKING/OBSOLETE)
# if (NOT DEFINED EXCLUDE_IAMCALC AND "iamcalc" IN_LIST EXECUTABLES)
#   add_executable_plus(iamcalc TRUE # This target uses BLAS
#       iamcalc.c
#       ftran.c
#       ulib.c
#       cnjj.c
#       slibgcc.c
#       catutil.c
#       lsqfit.c
#       # dblas.c removed
#   )
# endif()

# Add sortn
if (NOT DEFINED EXCLUDE_SORTN AND "sortn" IN_LIST EXECUTABLES)
  add_executable_plus(sortn FALSE # This target does NOT use BLAS
      sortn.c
      sortsub.c
  )
endif()

# Add reassign
if (NOT DEFINED EXCLUDE_REASSIGN AND "reassign" IN_LIST EXECUTABLES)
  add_executable_plus(reassign TRUE # This target uses BLAS
      reassign.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      # dblas.c removed
  )
endif()

# Add sortegy
if (NOT DEFINED EXCLUDE_SORTEGY AND "sortegy" IN_LIST EXECUTABLES)
  add_executable_plus(sortegy TRUE # This target uses BLAS
      sortegy.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      # dblas.c removed
  )
endif()

# Add iambak
if (NOT DEFINED EXCLUDE_IAMBAK AND "iambak" IN_LIST EXECUTABLES)
  add_executable_plus(iambak TRUE # This target uses BLAS
      iambak.c
      readopt.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      # dblas.c removed
  )
endif()