"""
Round-trip tests for TOML file format (Task 12 Phase 6).

mol.toml → load_fit_input → FitSession.from_input → run
         → save_fit_output → load_cat_input (+ int.toml) → CatSession → run
"""

import math
import pathlib
import tempfile
import pytest
import pickett
from pickett import (
    FitSession, CatSession,
    FitInput, CatInput,
    Parameter, LineRecord, DipoleMoment, VibState,
    SpinvOptions, DpiOptions, EngineOptions, EngineKind, CatControl,
    fit_input_to_dict, fit_input_from_dict,
    cat_input_to_dict,
    fit_output_to_dict, cat_output_to_dict,
    load_fit_input, load_cat_input,
    save_fit_output, save_cat_output,
    parse_fit_files, parse_cat_files,
)

REPO = pathlib.Path(__file__).parents[2]
FIXTURES = pathlib.Path(__file__).parent / "toml_fixtures"
CO4_BASE = str(REPO / "spfit_spcat_test_suite/diatomic_molecules/co_4/co_4")
CO4_VAR  = str(REPO / "spfit_spcat_test_suite/diatomic_molecules/co_4/v2008_results/co_4")


# ── Struct round-trip (no files) ────────────────────────────────────────────

class TestFitInputDictRoundTrip:
    def setup_method(self):
        self.fi = parse_fit_files(CO4_BASE + ".par", CO4_BASE + ".lin")
        # parse_fit_files uses raw_lines; copy them into lines for TOML round-trip
        fi2 = pickett.parse_fit_files(CO4_BASE + ".par", CO4_BASE + ".lin")
        # Build lines from raw_lines by running through the file-based path:
        # simplest: use the fixture file which has explicit lines
        self.fi_with_lines = load_fit_input(FIXTURES / "co_4.toml")

    def test_parameter_count(self):
        d = fit_input_to_dict(self.fi_with_lines)
        fi2 = fit_input_from_dict(d)
        assert len(fi2.parameters) == len(self.fi_with_lines.parameters)

    def test_parameter_values_preserved(self):
        d = fit_input_to_dict(self.fi_with_lines)
        fi2 = fit_input_from_dict(d)
        for p1, p2 in zip(self.fi_with_lines.parameters, fi2.parameters):
            assert p1.id == p2.id
            assert math.isclose(p1.value, p2.value, rel_tol=1e-14)
            assert math.isclose(p1.a_priori_error, p2.a_priori_error, rel_tol=1e-6)
            assert p1.fixed == p2.fixed
            assert p1.label == p2.label

    def test_engine_options_preserved(self):
        d = fit_input_to_dict(self.fi_with_lines)
        fi2 = fit_input_from_dict(d)
        e1 = self.fi_with_lines.engine_options
        e2 = fi2.engine_options
        assert e1.kind == e2.kind
        s1, s2 = e1.spinv, e2.spinv
        assert s1.nam_file == s2.nam_file
        assert len(s1.vibs) == len(s2.vibs)
        v1, v2 = s1.vibs[0], s2.vibs[0]
        assert v1.knmax == v2.knmax
        assert v1.esym_weight == v2.esym_weight
        assert v1.spin_degeneracies == v2.spin_degeneracies
        assert v1.symmetric_rotor_quanta == v2.symmetric_rotor_quanta

    def test_line_records_preserved(self):
        d = fit_input_to_dict(self.fi_with_lines)
        fi2 = fit_input_from_dict(d)
        assert len(fi2.lines) == len(self.fi_with_lines.lines)
        for lr1, lr2 in zip(self.fi_with_lines.lines, fi2.lines):
            assert lr1.nqn == lr2.nqn
            assert lr1.qn[:2*lr1.nqn] == lr2.qn[:2*lr2.nqn]
            assert math.isclose(lr1.freq, lr2.freq, rel_tol=1e-12)
            assert math.isclose(lr1.err, lr2.err, rel_tol=1e-6)

    def test_scalar_fields_preserved(self):
        d = fit_input_to_dict(self.fi_with_lines)
        fi2 = fit_input_from_dict(d)
        assert fi2.n_iterations == self.fi_with_lines.n_iterations
        assert fi2.max_lines == self.fi_with_lines.max_lines
        assert math.isclose(fi2.marquardt_param, self.fi_with_lines.marquardt_param, rel_tol=1e-6)


class TestCatInputDictRoundTrip:
    def setup_method(self):
        self.ci = parse_cat_files(CO4_VAR + ".var", CO4_BASE + ".int")

    def test_control_preserved(self):
        d = cat_input_to_dict(self.ci)
        # cat_input_to_dict only stores control + dipoles (no engine/params)
        ctrl_d = d.get("control", {})
        c1, c2 = self.ci.control, pickett.toml_io._cat_control_from_dict(ctrl_d)
        assert c1.itag == c2.itag
        assert math.isclose(c1.qrot, c2.qrot, rel_tol=1e-6)
        assert math.isclose(c1.fqmax, c2.fqmax, rel_tol=1e-6)
        assert math.isclose(c1.thrsh, c2.thrsh, rel_tol=1e-6)

    def test_dipoles_preserved(self):
        d = cat_input_to_dict(self.ci)
        dips = [pickett.toml_io._dipole_from_dict(dm) for dm in d.get("dipoles", [])]
        assert len(dips) == len(self.ci.dipoles)
        for dm1, dm2 in zip(self.ci.dipoles, dips):
            assert dm1.id == dm2.id
            assert math.isclose(dm1.value, dm2.value, rel_tol=1e-6)


# ── Fixture file tests ──────────────────────────────────────────────────────

class TestLoadFixtures:
    def test_load_fit_input_fixture(self):
        fi = load_fit_input(FIXTURES / "co_4.toml")
        assert fi.title == "CO; G. Winnewisser, et al."
        assert fi.n_iterations == 4
        assert fi.max_lines == 31
        assert len(fi.parameters) == 4
        assert len(fi.lines) == 31

    def test_fixture_B_parameter(self):
        fi = load_fit_input(FIXTURES / "co_4.toml")
        B = next(p for p in fi.parameters if p.id == 100)
        assert math.isclose(B.value, 57635.96804050475, rel_tol=1e-10)
        assert B.label == "B"

    def test_fixture_engine_options(self):
        fi = load_fit_input(FIXTURES / "co_4.toml")
        assert fi.engine_options.kind == EngineKind.Spinv
        spinv = fi.engine_options.spinv
        assert spinv.nam_file == "spinl.nam"
        assert len(spinv.vibs) == 1
        v = spinv.vibs[0]
        assert v.knmax == 0
        assert v.esym_weight == -1
        assert v.spin_degeneracies == [1]
        assert v.symmetric_rotor_quanta is True

    def test_fixture_first_line(self):
        fi = load_fit_input(FIXTURES / "co_4.toml")
        lr = fi.lines[0]
        assert lr.nqn == 1
        assert lr.qn[:2] == [1, 0]
        assert math.isclose(lr.freq, 115271.2018, rel_tol=1e-10)
        assert math.isclose(lr.err, 0.0005, rel_tol=1e-6)


# ── Fit round-trip via TOML fixture ────────────────────────────────────────

class TestFitViaTOML:
    def test_co4_fit_via_toml_matches_file_based(self):
        """TOML fixture → FitSession.from_input → run matches legacy file-based result."""
        ref = pickett.fit_files(CO4_BASE)
        fi = load_fit_input(FIXTURES / "co_4.toml")
        out = FitSession.from_input(fi).run()
        assert math.isclose(out.xsqbest, ref.xsqbest, rel_tol=1e-6)
        assert len(out.par) == 4
        for p_ref, p_out in zip(ref.par, out.par):
            assert math.isclose(p_ref, p_out, rel_tol=1e-10)

    def test_co4_fit_xsqbest_reasonable(self):
        fi = load_fit_input(FIXTURES / "co_4.toml")
        out = FitSession.from_input(fi).run()
        assert out.xsqbest < 1.0


# ── Full round-trip: fit → save fitted.toml → load → cat ────────────────────

class TestFullRoundTrip:
    def test_co4_full_round_trip(self, tmp_path):
        """TOML fit input → fit → save mol.fitted.toml → load → cat → check freq."""
        # Step 1: load FitInput from fixture and run fit
        fi = load_fit_input(FIXTURES / "co_4.toml")
        fit_out = FitSession.from_input(fi).run()

        # Step 2: save mol.fitted.toml
        var_toml = tmp_path / "co.fitted.toml"
        save_fit_output(fit_out, fi, var_toml)
        assert var_toml.exists()

        # Step 3: load CatInput from fitted.toml + dipoles.toml fixture
        int_toml = FIXTURES / "co_4.dipoles.toml"
        ci = load_cat_input(var_toml, int_toml)
        assert len(ci.parameters) == 4
        assert ci.control.itag == 28503

        # Step 4: run catalog generation
        cat_out = CatSession.from_input(ci).run()
        assert cat_out.nline > 0

        # Step 5: check first catalog line frequency matches reference
        freq_str = cat_out.cat_lines[0][:13].strip()
        freq = float(freq_str)
        assert math.isclose(freq, 115271.2021, rel_tol=1e-4)

    def test_co4_var_toml_fields(self, tmp_path):
        """Saved mol.fitted.toml has expected top-level fields."""
        import tomllib as tl  # stdlib 3.11+; falls back via conftest

        fi = load_fit_input(FIXTURES / "co_4.toml")
        fit_out = FitSession.from_input(fi).run()
        var_toml = tmp_path / "co.fitted.toml"
        save_fit_output(fit_out, fi, var_toml)

        with open(var_toml, "rb") as f:
            try:
                d = tl.load(f)
            except Exception:
                import tomli as tl2  # type: ignore[import]
                f.seek(0)
                d = tl2.load(f)

        assert "title" in d
        assert "xsqbest" in d
        assert "parameters" in d
        assert "engine_options" in d
        assert len(d["parameters"]) == 4
        # check fitted B value is close to expected
        B_par = next(p for p in d["parameters"] if p["id"] == 100)
        assert math.isclose(B_par["value"], 57635.968040505, rel_tol=1e-8)

    def test_cat_output_saved(self, tmp_path):
        """save_cat_output produces a loadable TOML file."""
        import tomllib as tl
        ci = parse_cat_files(CO4_VAR + ".var", CO4_BASE + ".int")
        cat_out = CatSession.from_input(ci).run()
        cat_toml = tmp_path / "co.catalog.toml"
        save_cat_output(cat_out, cat_toml)

        with open(cat_toml, "rb") as f:
            try:
                d = tl.load(f)
            except Exception:
                import tomli as tl2  # type: ignore[import]
                f.seek(0)
                d = tl2.load(f)

        assert d["nline"] > 0
        assert len(d["cat_lines"]) == d["nline"]
        assert len(d["temp"]) > 0
        # first line frequency
        freq = float(d["cat_lines"][0][:13].strip())
        assert math.isclose(freq, 115271.2021, rel_tol=1e-4)
