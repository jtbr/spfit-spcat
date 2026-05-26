"""
File-free path tests for o2_1and3 (open-shell linear, triplet) and fso3
(C3v symmetric top, two vibronic states).

Each test uses parse_fit_files / parse_cat_files to build a FitInput / CatInput
struct, then drives FitSession.from_input() / CatSession.from_input() and
checks that the output is bit-identical to the file-based path and consistent
with the v2008 reference baseline.
"""

import math
import pickett
from pickett import (
    FitSession, CatSession,
    parse_fit_files, parse_cat_files,
    fit_files, cat_from_files,
)

O2_BASE = "spfit_spcat_test_suite/diatomic_molecules/o2_1and3/o2_1and3"
O2_VAR  = "spfit_spcat_test_suite/diatomic_molecules/o2_1and3/v2008_results/o2_1and3"

FSO3_BASE = "spfit_spcat_test_suite/symmetric_tops/fso3/fso3"
FSO3_VAR  = "spfit_spcat_test_suite/symmetric_tops/fso3/v2008_results/fso3"


# ── Regression test: static-cache bug in setgsym ─────────────────────────────
#
# setgsym() used to cache the last gsym in a function-local static.  When two
# CalFit/CalCat instances ran sequentially in the same process with the same
# gsym value, the second instance's SpinvContext was never initialized —
# causing completely wrong energy levels and "NEXT LINE NOT USED IN FIT" for
# most transitions.  The fix removes the static so every ctx is initialized.

def test_o2_fit_files_second_run_gives_same_result():
    """fit_files for o2 twice in the same process must give identical xsqbest."""
    r1 = fit_files(O2_BASE)
    r2 = fit_files(O2_BASE)
    assert abs(r1.xsqbest - r2.xsqbest) < 1e-8, (
        f"Second fit_files call gave different xsqbest: {r1.xsqbest} vs {r2.xsqbest}"
    )

def test_fso3_fit_files_second_run_gives_same_result():
    """fit_files for fso3 twice in the same process must give identical xsqbest."""
    r1 = fit_files(FSO3_BASE)
    r2 = fit_files(FSO3_BASE)
    assert abs(r1.xsqbest - r2.xsqbest) < 1e-8, (
        f"Second fit_files call gave different xsqbest: {r1.xsqbest} vs {r2.xsqbest}"
    )


# ── O2 — fit ─────────────────────────────────────────────────────────────────

class TestO2Fit:
    def test_round_trip_matches_file_based(self):
        """parse_fit_files → from_input → run is bit-identical to fit_files."""
        fi  = parse_fit_files(O2_BASE + ".par", O2_BASE + ".lin")
        out = FitSession.from_input(fi).run()
        ref = fit_files(O2_BASE)

        assert abs(ref.xsqbest - out.xsqbest) < 1e-8
        assert ref.par   == out.par
        assert ref.erpar == out.erpar
        assert ref.itr   == out.itr

    def test_xsqbest_matches_reference(self):
        """xsqbest ≈ 0.88877 (v2008 reference, 1 iteration)."""
        fi  = parse_fit_files(O2_BASE + ".par", O2_BASE + ".lin")
        out = FitSession.from_input(fi).run()

        assert math.isfinite(out.xsqbest)
        assert abs(out.xsqbest - 0.88877) < 1e-4
        assert out.itr == 1

    def test_parameter_count_and_B(self):
        """O2 has 13 parameters; B ≈ 43100.44276 MHz."""
        fi  = parse_fit_files(O2_BASE + ".par", O2_BASE + ".lin")
        out = FitSession.from_input(fi).run()

        assert len(out.par) == 13
        # B is parameter index 1 (index 0 is E(0))
        B_ref = 4.310044276e4   # from v2008 .fit
        assert math.isclose(out.par[1], B_ref, rel_tol=1e-6)


# ── O2 — cat ─────────────────────────────────────────────────────────────────

class TestO2Cat:
    def test_round_trip_matches_file_based(self):
        """parse_cat_files → from_input → run matches cat_from_files output."""
        ci  = parse_cat_files(O2_VAR + ".var", O2_BASE + ".int")
        out = CatSession.from_input(ci).run()
        ref = cat_from_files(O2_BASE + ".int", O2_VAR + ".var")

        assert len(out.cat_lines) == len(ref.cat_lines)
        assert out.cat_lines[:5]  == ref.cat_lines[:5]
        assert out.qsum           == ref.qsum

    def test_catalog_nonempty_and_sorted(self):
        """O2 catalog has lines and they are sorted by frequency."""
        ci  = parse_cat_files(O2_VAR + ".var", O2_BASE + ".int")
        out = CatSession.from_input(ci).run()

        assert out.nline > 0
        freqs = [float(line[:13]) for line in out.cat_lines]
        assert freqs == sorted(freqs)


# ── FSO3 — fit ───────────────────────────────────────────────────────────────

class TestFso3Fit:
    def test_round_trip_matches_file_based(self):
        """parse_fit_files → from_input → run is bit-identical to fit_files."""
        fi  = parse_fit_files(FSO3_BASE + ".par", FSO3_BASE + ".lin")
        out = FitSession.from_input(fi).run()
        ref = fit_files(FSO3_BASE)

        assert abs(ref.xsqbest - out.xsqbest) < 1e-8
        assert ref.par   == out.par
        assert ref.erpar == out.erpar
        assert ref.itr   == out.itr

    def test_xsqbest_matches_reference(self):
        """xsqbest ≈ 0.93610 (v2008 reference, 1 iteration)."""
        fi  = parse_fit_files(FSO3_BASE + ".par", FSO3_BASE + ".lin")
        out = FitSession.from_input(fi).run()

        assert math.isfinite(out.xsqbest)
        assert abs(out.xsqbest - 0.93610) < 1e-4
        assert out.itr == 1

    def test_parameter_count_and_B(self):
        """FSO3 has 17 parameters; B ≈ 5195.528 MHz."""
        fi  = parse_fit_files(FSO3_BASE + ".par", FSO3_BASE + ".lin")
        out = FitSession.from_input(fi).run()

        assert len(out.par) == 17
        # B (parameter id 199) is index 1 in the .par file
        B_ref = 5.195528437873451e3   # from v2008 .fit
        assert math.isclose(out.par[1], B_ref, rel_tol=1e-6)


# ── FSO3 — cat ───────────────────────────────────────────────────────────────

class TestFso3Cat:
    def test_round_trip_matches_file_based(self):
        """parse_cat_files → from_input → run matches cat_from_files output."""
        ci  = parse_cat_files(FSO3_VAR + ".var", FSO3_BASE + ".int")
        out = CatSession.from_input(ci).run()
        ref = cat_from_files(FSO3_BASE + ".int", FSO3_VAR + ".var")

        assert len(out.cat_lines) == len(ref.cat_lines)
        assert out.cat_lines[:5]  == ref.cat_lines[:5]
        assert out.qsum           == ref.qsum

    def test_catalog_nonempty_and_sorted(self):
        """FSO3 catalog has lines and they are sorted by frequency."""
        ci  = parse_cat_files(FSO3_VAR + ".var", FSO3_BASE + ".int")
        out = CatSession.from_input(ci).run()

        assert out.nline > 0
        freqs = [float(line[:13]) for line in out.cat_lines]
        assert freqs == sorted(freqs)
