cmake_minimum_required(VERSION 3.0)
project(spfit)
cmake_policy(SET CMP0057 NEW)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -march=native")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O3 -march=native")

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

# Find BLAS
find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    include_directories(${BLAS_INCLUDE_DIRS})
    set(BLAS_LIBRARIES ${BLAS_LIBRARIES})
else()
    message(STATUS "BLAS not found, using dblas fallback")
    set(BLAS_LIBRARIES dblas)
endif()

# Define source files
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
    calcat.c
    sortsub.c
    calmrg.c
    calbak.c
    termval.c
    stark.c
    moiam.c
    iamcalc.c
    cnvwn.c
    sortn.c
    reassign.c
    sortegy.c
    iambak.c
    ulib.c
    cnjj.c
    slibgcc.c
    lsqfit.c
    dblas.c
    dpi.c
    ftran.c
)

# Add executable macro
macro(add_executable_plus name sources)
  add_executable(${name} ${sources})
  if ("${name}" MATCHES "^(spfit|spcat)$")
    target_link_libraries(${name} ${BLAS_LIBRARIES} -lm -lstdc++)
  else()
    target_link_libraries(${name} ${BLAS_LIBRARIES} -lm)
  endif()
  set_property(TARGET ${name} PROPERTY LINKER_LANGUAGE CXX)
endmacro()

# Define executables
set(DEFAULT_EXECUTABLES spfit spcat)
set(EXECUTABLES ${DEFAULT_EXECUTABLES})

set_source_files_properties(${SOURCE_FILES} PROPERTIES COMPILE_FLAGS "-Wno-unused-function")

project(
  Pickett
  DESCRIPTION "SPFIT/SPCAT"
  LANGUAGES CXX C)

# Add spfit
add_executable_plus(spfit
    calfit.cpp
    subfit.c
    spinv_setup.c
    spinv_spin_symmetry.c
    spinv_linalg_sort.c
    spinv_hamiltonian.c
    spinv_utils.c
    spinit.c
    SpinvEngine.cpp
    catutil.c
    ulib.c
    cnjj.c
    slibgcc.c
    lsqfit.c
    dblas.c
    dpi.c
    ftran.c
)

# Add spcat
add_executable_plus(spcat
    calcat.c
    sortsub.c
    spinv_setup.c
    spinv_spin_symmetry.c
    spinv_linalg_sort.c
    spinv_hamiltonian.c
    spinv_utils.c
    spinit.c
    SpinvEngine.cpp
    catutil.c
    ulib.c
    cnjj.c
    slibgcc.c
    lsqfit.c
    dblas.c
)

# Add calmrg
if (NOT DEFINED EXCLUDE_CALMRG AND "calmrg" IN_LIST [${EXECUTABLES}])
  add_executable_plus(calmrg
      calmrg.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      dblas.c
  )
endif()

# Add calbak
if (NOT DEFINED EXCLUDE_CALBAK AND "calbak" IN_LIST [${EXECUTABLES}])
  add_executable_plus(calbak
      calbak.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      dblas.c
  )
endif()

# Add termval
if (NOT DEFINED EXCLUDE_TERMVAL AND "termval" IN_LIST [${EXECUTABLES}])
  add_executable_plus(termval
      termval.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      dblas.c
  )
endif()

# Add stark
if (NOT DEFINED EXCLUDE_STARK AND "stark" IN_LIST [${EXECUTABLES}])
  add_executable_plus(stark
      stark.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      dblas.c
  )
endif()

# Add moiam
if (NOT DEFINED EXCLUDE_MOIAM AND "moiam" IN_LIST [${EXECUTABLES}])
  add_executable_plus(moiam
      moiam.c
      ftran.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      dblas.c
  )
endif()

# # Add iamcalc (NOT WORKING/OBSOLETE)
# if (NOT DEFINED EXCLUDE_IAMCALC AND iamcalc IN_LIST [${EXECUTABLES}])
#   add_executable_plus(iamcalc
#       iamcalc.c
#       ftran.c
#       ulib.c
#       cnjj.c
#       slibgcc.c
#       catutil.c
#       lsqfit.c
#       dblas.c
#   )
# endif()

# Add sortn
if (NOT DEFINED EXCLUDE_SORTN AND sortn IN_LIST [${EXECUTABLES}])
  add_executable_plus(sortn
      sortn.c
      sortsub.c
  )
endif()

# Add reassign
if (NOT DEFINED EXCLUDE_REASSIGN AND reassign IN_LIST [${EXECUTABLES}])
  add_executable_plus(reassign
      reassign.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      dblas.c
  )
endif()

# Add sortegy
if (NOT DEFINED EXCLUDE_SORTEGY AND sortegy IN_LIST [${EXECUTABLES}])
  add_executable_plus(sortegy
      sortegy.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      dblas.c
  )
endif()

# Add iambak
if (NOT DEFINED EXCLUDE_IAMBAK AND iambak IN_LIST [${EXECUTABLES}])
  add_executable_plus(iambak
      iambak.c
      readopt.c
      ulib.c
      cnjj.c
      slibgcc.c
      catutil.c
      lsqfit.c
      dblas.c
  )
endif()
