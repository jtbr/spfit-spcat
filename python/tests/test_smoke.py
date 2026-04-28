"""
Smoke tests: round-trip co_4 through the Python API and compare key output
values against the reference files in spfit_spcat_test_suite/.
"""
import pathlib
import math
import pytest
import pickett

# Locate the co_4 test molecule relative to this file's position in the repo.
REPO_ROOT = pathlib.Path(__file__).parents[2]
CO4_BASE  = str(REPO_ROOT / "spfit_spcat_test_suite/diatomic_molecules/co_4/co_4")
CO4_REF   = REPO_ROOT / "spfit_spcat_test_suite/diatomic_molecules/co_4/reference_outputs"
# The .var file needed by cat_files comes from the pre-computed v2008 results.
CO4_VAR_BASE = str(REPO_ROOT / "spfit_spcat_test_suite/diatomic_molecules/co_4/v2008_results/co_4")


# ---------------------------------------------------------------------------
# CalFit (high-level)
# ---------------------------------------------------------------------------

class TestFitFiles:
    def setup_method(self):
        self.out = pickett.fit_files(CO4_BASE)

    def test_returns_calfit_output(self):
        assert isinstance(self.out, pickett.CalFitOutput)

    def test_parameter_count(self):
        # co_4 has 4 fitted parameters
        assert len(self.out.par) == 4
        assert len(self.out.erpar) == 4

    def test_xsqbest_reasonable(self):
        # Reference fit achieves near-zero rms for co_4
        assert self.out.xsqbest < 1.0

    def test_B_constant(self):
        # First parameter is B ≈ 57635.968 MHz; allow 1 ppm tolerance
        B_ref = 5.7635968040505e4  # from reference_outputs/co_4.fit
        assert math.isclose(self.out.par[0], B_ref, rel_tol=1e-6)

    def test_iterations(self):
        assert self.out.itr >= 1


# ---------------------------------------------------------------------------
# CalFit (low-level session path)
# ---------------------------------------------------------------------------

class TestFitSession:
    def test_session_round_trip(self):
        session = pickett.FitSession(CO4_BASE + ".par", CO4_BASE + ".lin")
        out = session.run()
        assert isinstance(out, pickett.CalFitOutput)
        assert len(out.par) == 4
        assert out.xsqbest < 1.0

    def test_session_single_use(self):
        session = pickett.FitSession(CO4_BASE + ".par", CO4_BASE + ".lin")
        session.run()
        with pytest.raises(pickett.ValidationError):
            session.run()

    def test_bad_engine_raises(self):
        with pytest.raises(pickett.ValidationError):
            pickett.FitSession(CO4_BASE + ".par", CO4_BASE + ".lin", engine="bad")


# ---------------------------------------------------------------------------
# CalCat (high-level)
# ---------------------------------------------------------------------------

class TestCatFiles:
    def setup_method(self):
        # int file from co_4, var from pre-computed v2008 result
        self.out = pickett.cat_files(CO4_VAR_BASE, int_path=CO4_BASE)

    def test_returns_calcat_output(self):
        assert isinstance(self.out, pickett.CalCatOutput)

    def test_catalog_lines_nonempty(self):
        assert self.out.nline > 0
        assert len(self.out.cat_lines) > 0

    def test_first_catalog_line_frequency(self):
        # co_4 reference: first line near 115271.2021 MHz
        # Catalog format: cols 0-12 are frequency in MHz (right-justified)
        freq_str = self.out.cat_lines[0][:13].strip()
        freq = float(freq_str)
        assert math.isclose(freq, 115271.2021, rel_tol=1e-5)

    def test_catalog_sorted_by_frequency(self):
        freqs = [float(line[:13]) for line in self.out.cat_lines]
        assert freqs == sorted(freqs)

    def test_partition_function_300K(self):
        # Reference: Q(300 K) should be a reasonable positive number for CO
        temps = self.out.temp
        qsums = self.out.qsum
        q300 = None
        best = 1e30
        for t, q in zip(temps, qsums):
            if abs(t - 300.0) < best:
                best = abs(t - 300.0)
                q300 = q
        assert q300 is not None and q300 > 0

    def test_partition_function_table_length(self):
        assert self.out.ntemp > 0
        assert len(self.out.temp) == self.out.ntemp
        assert len(self.out.qsum) == self.out.ntemp


# ---------------------------------------------------------------------------
# CalCat (low-level session path)
# ---------------------------------------------------------------------------

class TestCatSession:
    def test_session_round_trip(self):
        session = pickett.CatSession(CO4_BASE + ".int", CO4_VAR_BASE + ".var")
        out = session.run()
        assert isinstance(out, pickett.CalCatOutput)
        assert out.nline > 0

    def test_session_single_use(self):
        session = pickett.CatSession(CO4_BASE + ".int", CO4_VAR_BASE + ".var")
        session.run()
        with pytest.raises(pickett.ValidationError):
            session.run()


# ---------------------------------------------------------------------------
# Exception mapping
# ---------------------------------------------------------------------------

class TestExceptions:
    def test_io_error_on_missing_file(self):
        with pytest.raises(pickett.IoError):
            pickett.fit_files("/nonexistent/path/molecule")

    def test_io_error_is_cal_error(self):
        with pytest.raises(pickett.CalError):
            pickett.fit_files("/nonexistent/path/molecule")
