"""
Tests for typed-struct input API:
  parse_fit_files / parse_cat_files → FitInput / CatInput
  FitSession.from_input / CatSession.from_input
"""

import math
import pytest
import pickett
from pickett import (
    FitSession, CatSession, FitInput, CatInput,
    Parameter, LineRecord, DipoleMoment,
    SpinvOptions, EngineOptions, EngineKind,
    parse_fit_files, parse_cat_files,
    fit_files, cat_from_files,
    InputError,
)

CO4_BASE = "spfit_spcat_test_suite/diatomic_molecules/co_4/co_4"
CO4_VAR  = "spfit_spcat_test_suite/diatomic_molecules/co_4/v2008_results/co_4"


# ── parse_fit_files round-trip ────────────────────────────────────────────────

def test_fit_round_trip_matches_file_based():
    """parse_fit_files → FitSession.from_input → run matches fit_files result."""
    ref = fit_files(CO4_BASE)
    fi = parse_fit_files(CO4_BASE + ".par", CO4_BASE + ".lin")
    out = FitSession.from_input(fi).run()

    assert abs(ref.xsqbest - out.xsqbest) < 1e-8
    assert ref.par == out.par
    assert ref.erpar == out.erpar
    assert ref.itr == out.itr


def test_parse_fit_files_struct_contents():
    """parse_fit_files returns well-formed FitInput for CO v=0."""
    fi = parse_fit_files(CO4_BASE + ".par", CO4_BASE + ".lin")

    # CO has 4 parameters: B, D, H, L
    assert len(fi.parameters) == 4
    # B rotational constant ≈ 57635.968 MHz
    B = next(p for p in fi.parameters if p.id == 100)
    assert abs(B.value - 57635.968040505) < 1e-3
    assert B.a_priori_error > 1e30   # large (free to fit)
    assert not B.fixed

    # lines list is empty — parse_fit_files uses an internal raw-lines fast path
    assert fi.lines == []


# ── parse_cat_files round-trip ────────────────────────────────────────────────

def test_cat_round_trip_matches_file_based():
    """parse_cat_files → CatSession.from_input → run matches cat_from_files result."""
    ref = cat_from_files(CO4_BASE + ".int", CO4_VAR + ".var")
    ci = parse_cat_files(CO4_VAR + ".var", CO4_BASE + ".int")
    out = CatSession.from_input(ci).run()

    assert len(out.cat_lines) == len(ref.cat_lines)
    assert out.cat_lines[:5] == ref.cat_lines[:5]
    assert out.qsum == ref.qsum


# ── manual FitInput construction ──────────────────────────────────────────────

def test_fit_modify_parameter_then_run():
    """Modifying a parsed FitInput's parameter value and re-running works."""
    fi = parse_fit_files(CO4_BASE + ".par", CO4_BASE + ".lin")

    # Zero out the D (centrifugal distortion) parameter so it doesn't fit
    for p in fi.parameters:
        if p.id == 200:   # -D
            p.value = 0.0
            p.a_priori_error = -1.0   # fix it (negative → dependent)
            break

    out = FitSession.from_input(fi).run()
    # Result should differ from default (3-param fit, fewer dof)
    # but should complete without error and give a finite xsqbest
    assert math.isfinite(out.xsqbest)
    assert len(out.par) == 4


# ── Parameter struct field validation ────────────────────────────────────────

def test_parameter_fields():
    fi = parse_fit_files(CO4_BASE + ".par", CO4_BASE + ".lin")
    p = fi.parameters[0]  # B

    assert isinstance(p.id, int)
    assert isinstance(p.value, float)
    assert isinstance(p.a_priori_error, float)
    assert isinstance(p.fixed, bool)
    assert isinstance(p.label, str)


def test_fit_input_empty_lines_raises():
    fi = FitInput()
    fi.parameters = [Parameter()]
    # fi.lines stays empty; no raw_lines (internal field) populated
    with pytest.raises(InputError):
        FitSession.from_input(fi)


# ── EngineKind enum ───────────────────────────────────────────────────────────

def test_engine_kind_values():
    assert EngineKind.Spinv != EngineKind.Dpi


def test_parse_fit_files_default_engine_is_spinv():
    fi = parse_fit_files(CO4_BASE + ".par", CO4_BASE + ".lin")
    assert fi.engine_options.kind == EngineKind.Spinv
