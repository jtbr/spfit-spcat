cmake_minimum_required(VERSION 3.15)
project(pickett_python LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)

# -ffp-contract=off and bundled dblas.c are required for numerical reproducibility
# (see project CLAUDE.md for details on FMA / Cholesky sign-flip issue).
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffp-contract=off -Wno-missing-field-initializers")
set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3 -march=native -ffp-contract=off -Wno-missing-field-initializers")

find_package(Python 3.8 COMPONENTS Interpreter Development.Module REQUIRED)

# Locate nanobind cmake files via the Python interpreter so the same venv is used
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
find_package(nanobind CONFIG REQUIRED PATHS "${NB_DIR}")

set(PARENT_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../src")
include_directories("${PARENT_SRC}")

set(ENGINE_SRCS
    ${PARENT_SRC}/engine/spinv_setup.cpp
    ${PARENT_SRC}/engine/spinv_spin_symmetry.cpp
    ${PARENT_SRC}/engine/spinv_linalg_sort.c
    ${PARENT_SRC}/engine/spinv_hamiltonian.cpp
    ${PARENT_SRC}/engine/spinv_utils.cpp
    ${PARENT_SRC}/engine/spinit.cpp
    ${PARENT_SRC}/engine/dpi.cpp
    ${PARENT_SRC}/engine/SpinvEngine.cpp
    ${PARENT_SRC}/engine/DpiEngine.cpp
)

set(SPLIB_SRCS
    ${PARENT_SRC}/splib/catutil.c
    ${PARENT_SRC}/splib/ulib.c
    ${PARENT_SRC}/splib/cnjj.c
    ${PARENT_SRC}/splib/lsqfit.c
    ${PARENT_SRC}/splib/ftran.c
    ${PARENT_SRC}/splib/dblas.c
)

set(COMMON_SRCS
    ${PARENT_SRC}/common/Logger.cpp
    ${PARENT_SRC}/common/SigintFlag.cpp
    ${PARENT_SRC}/common/file_helpers.cpp
)

set(SPFIT_SRCS
    ${PARENT_SRC}/spfit/CalFit.cpp
    ${PARENT_SRC}/spfit/CalFit_helpers.cpp
    ${PARENT_SRC}/spfit/CalFitIO.cpp
    ${PARENT_SRC}/spfit/subfit.cpp
)

set(SPCAT_SRCS
    ${PARENT_SRC}/spcat/CalCat.cpp
    ${PARENT_SRC}/spcat/CalCat_helpers.cpp
    ${PARENT_SRC}/spcat/CalCatIO.cpp
)

nanobind_add_module(_pickett NB_STATIC
    src/bindings.cpp
    ${ENGINE_SRCS}
    ${SPLIB_SRCS}
    ${COMMON_SRCS}
    ${SPFIT_SRCS}
    ${SPCAT_SRCS}
)

target_link_libraries(_pickett PRIVATE m)
install(TARGETS _pickett LIBRARY DESTINATION pickett)
