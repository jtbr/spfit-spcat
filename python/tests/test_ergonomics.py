"""
Tests for the ergonomic surface of the typed-struct API:
  - keyword-argument constructors for every input struct
  - LineRecord.qn pad-on-write / Sequence acceptance
  - CatControl.iflg convenience properties
  - EngineOptions auto-inference of `kind` from spinv/dpi args
  - __repr__ output is informative
"""

import pytest
import pickett
from pickett import (
    Parameter, LineRecord, DipoleMoment, VibState,
    SpinvOptions, DpiOptions, EngineOptions, EngineKind,
    CatControl, FitInput, CatInput, FitSession,
)


# ── Keyword constructors ──────────────────────────────────────────────────────

def test_parameter_kwargs():
    p = Parameter(id=100, value=57635.968, a_priori_error=1e37, fixed=False, label="B")
    assert p.id == 100
    assert p.value == 57635.968
    assert p.a_priori_error == 1e37
    assert p.fixed is False
    assert p.label == "B"


def test_parameter_defaults():
    p = Parameter()
    assert p.id == 0 and p.value == 0.0 and p.a_priori_error == 1e37
    assert p.fixed is False and p.label == ""


def test_line_record_kwargs():
    lr = LineRecord(qn=[1, 0], nqn=1, freq=115271.2018, err=0.05, weight=1.0)
    assert lr.qn[0] == 1 and lr.qn[1] == 0 and lr.qn[2] == 0  # padded
    assert lr.nqn == 1
    assert lr.freq == 115271.2018
    assert lr.err == 0.05


def test_dipole_kwargs():
    d = DipoleMoment(id=11, value=0.112, starts_new_component=True)
    assert d.id == 11 and d.value == 0.112 and d.starts_new_component is True


def test_vib_state_kwargs():
    v = VibState(knmin=0, knmax=0, symmetric_rotor_quanta=True,
                 spin_degeneracies=[3], esym_weight=99)
    assert v.knmin == 0 and v.knmax == 0
    assert v.symmetric_rotor_quanta is True
    assert v.spin_degeneracies == [3]
    assert v.esym_weight == 99


def test_spinv_options_kwargs():
    so = SpinvOptions(oblate=True, vibs=[VibState(), VibState(index=1)])
    assert so.oblate is True
    assert len(so.vibs) == 2 and so.vibs[1].index == 1


def test_engine_options_auto_kind_from_spinv():
    eo = EngineOptions(spinv=SpinvOptions(vibs=[VibState()]))
    assert eo.kind == EngineKind.Spinv


def test_engine_options_auto_kind_from_dpi():
    eo = EngineOptions(dpi=DpiOptions(isdgn=2, nvib=3))
    assert eo.kind == EngineKind.Dpi
    assert eo.dpi.isdgn == 2 and eo.dpi.nvib == 3


def test_engine_options_default_kind_is_spinv():
    eo = EngineOptions()
    assert eo.kind == EngineKind.Spinv


def test_engine_options_explicit_kind_overrides():
    eo = EngineOptions(kind=EngineKind.Dpi, spinv=SpinvOptions())
    assert eo.kind == EngineKind.Dpi   # explicit kind wins


def test_cat_control_kwargs():
    cc = CatControl(itag=28503, fqmax=500.0, tmq=300.0, maxv=2)
    assert cc.itag == 28503 and cc.fqmax == 500.0 and cc.maxv == 2


def test_fit_input_kwargs():
    fi = FitInput(
        title="t",
        n_iterations=3,
        parameters=[Parameter(id=100, value=1.0)],
        lines=[LineRecord(qn=[1, 0], nqn=1, freq=100.0)],
    )
    assert fi.title == "t" and fi.n_iterations == 3
    assert len(fi.parameters) == 1 and len(fi.lines) == 1


def test_cat_input_kwargs():
    ci = CatInput(
        title="cat",
        control=CatControl(itag=42),
        dipoles=[DipoleMoment(id=1, value=0.5)],
    )
    assert ci.title == "cat" and ci.control.itag == 42
    assert len(ci.dipoles) == 1


# ── LineRecord.qn pad / accept Sequence ───────────────────────────────────────

def test_qn_short_list_pads_with_zeros():
    lr = LineRecord()
    lr.qn = [1, 0]
    assert lr.qn[:5] == [1, 0, 0, 0, 0]
    assert len(lr.qn) == 20
    assert all(v == 0 for v in lr.qn[2:])


def test_qn_empty_list_is_zero_padded():
    lr = LineRecord()
    lr.qn = []
    assert lr.qn == [0] * 20


def test_qn_too_long_raises():
    lr = LineRecord()
    with pytest.raises(ValueError, match="qn list too long"):
        lr.qn = [0] * 21


def test_qn_tuple_accepted():
    lr = LineRecord()
    lr.qn = (1, 0, 2, 3)
    assert lr.qn[:5] == [1, 0, 2, 3, 0]


# ── CatControl iflg convenience properties ────────────────────────────────────

def test_cc_output_strengths_sets_flag():
    cc = CatControl()
    cc.output_strengths = 1
    assert cc.iflg == 10
    cc.output_strengths = 2
    assert cc.iflg == 20


def test_cc_output_energies_sets_flag():
    cc = CatControl()
    cc.output_energies = 3
    assert cc.iflg == 3


def test_cc_wavenumbers_sets_flag():
    cc = CatControl()
    cc.wavenumbers = True
    assert cc.iflg == 1000
    cc.wavenumbers = False
    assert cc.iflg == 0


def test_cc_flags_compose():
    cc = CatControl()
    cc.wavenumbers = True
    cc.output_strengths = 1
    cc.output_energies = 2
    assert cc.iflg == 1012


def test_cc_flags_round_trip():
    cc = CatControl(iflg=1023)
    assert cc.wavenumbers is True
    assert cc.output_strengths == 2
    assert cc.output_energies == 3


# ── __repr__ smoke checks (just ensure they don't crash and contain key info) ─

def test_repr_parameter():
    p = Parameter(id=100, value=57635.968040505, label="B")
    s = repr(p)
    assert "Parameter(" in s
    assert "id=100" in s
    assert "57635.968" in s    # precision preserved
    assert "label='B'" in s


def test_repr_line_record_shows_only_meaningful_qn():
    lr = LineRecord(qn=[1, 0], nqn=1, freq=115271.2)
    s = repr(lr)
    # only the 2*nqn meaningful slots should appear
    assert "qn=[1, 0]" in s
    assert "qn=[1, 0, 0, 0" not in s


def test_repr_engine_options_includes_kind():
    eo = EngineOptions(dpi=DpiOptions())
    assert "EngineKind.Dpi" in repr(eo)


# ── End-to-end: kwarg construction → fit ──────────────────────────────────────

def test_full_co_fit_with_kwargs():
    """The quick-start example from API.md must run end-to-end."""
    fi = FitInput(
        title="CO v=0",
        n_iterations=5,
        parameters=[
            Parameter(id=100, value=57635.968, a_priori_error=1e37, label="B"),
            Parameter(id=200, value=-0.184,    a_priori_error=1e37, label="-D"),
            Parameter(id=300, value=1.7e-9,    a_priori_error=1e37, label="H"),
            Parameter(id=400, value=0.0,       fixed=True,           label="L"),
        ],
        lines=[
            LineRecord(qn=[1, 0], nqn=1, freq=115271.2018, err=0.05),
            LineRecord(qn=[2, 1], nqn=1, freq=230538.0000, err=0.05),
            LineRecord(qn=[3, 2], nqn=1, freq=345795.9899, err=0.05),
            LineRecord(qn=[4, 3], nqn=1, freq=461040.7681, err=0.05),
        ],
        engine_options=EngineOptions(spinv=SpinvOptions(vibs=[
            VibState(knmin=0, knmax=0, symmetric_rotor_quanta=True),
        ])),
    )
    out = FitSession.from_input(fi).run()
    assert out.itr >= 1
    assert len(out.par) == 4
    # B should be near its starting value (CO is well-determined)
    assert abs(out.par[0] - 57635.968) < 1.0
