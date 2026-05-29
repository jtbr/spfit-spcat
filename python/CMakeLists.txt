cmake_minimum_required(VERSION 3.15)
project(pickett_python LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_STANDARD 11)

# -ffp-contract=off / /fp:precise and bundled dblas.c are required for numerical reproducibility
# (see project CLAUDE.md for details on FMA / Cholesky sign-flip issue).
# Note: -march=native is intentionally absent here so wheels run on any compatible CPU.
if(MSVC)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)
    string(APPEND CMAKE_CXX_FLAGS " /O2 /fp:precise")
    string(APPEND CMAKE_C_FLAGS   " /O2 /fp:precise")
else()
    string(APPEND CMAKE_CXX_FLAGS " -O3 -ffp-contract=off -Wno-missing-field-initializers")
    string(APPEND CMAKE_C_FLAGS   " -O3 -ffp-contract=off -Wno-missing-field-initializers")
endif()

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

set(API_SRCS
    ${PARENT_SRC}/api/builders.cpp
    ${PARENT_SRC}/api/legacy_parser.cpp
)

nanobind_add_module(_pickett NB_STATIC
    src/bindings.cpp
    ${ENGINE_SRCS}
    ${SPLIB_SRCS}
    ${COMMON_SRCS}
    ${SPFIT_SRCS}
    ${SPCAT_SRCS}
    ${API_SRCS}
)

if(NOT MSVC)
    target_link_libraries(_pickett PRIVATE m)
endif()
install(TARGETS _pickett LIBRARY DESTINATION pickett)
