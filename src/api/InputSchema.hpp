/**
 * InputSchema.hpp — typed input records for CalFit / CalCat without filesystem I/O.
 *
 * These structs are the canonical public API.  The legacy file-reading path
 * (CalFitIO, CalCatIO) parses .par/.lin/.int/.var text into these structs and
 * then calls build_fit_input / build_cat_input (src/api/builders.hpp).
 *
 * Field naming follows the SPFIT/SPCAT .par/.int file specifications documented
 * in spinv.md and dpi.md.
 */

#ifndef INPUT_SCHEMA_HPP
#define INPUT_SCHEMA_HPP

#include <array>
#include <cstdint>
#include <string>
#include <vector>
#include "splib/calpgm_types.h"  // for MAXQN, MAXCAT

// ── Parameter ──────────────────────────────────────────────────────────────

struct Parameter {
    int64_t id = 0;                // decimal parameter ID; builder BCD-packs it
    double  value = 0.0;           // absolute value in MHz or cm-1
    double  a_priori_error = 1e37; // 1-sigma uncertainty; <= 1e-37 effectively fixes parameter
    bool    fixed = false;         // if true, excluded from fit (sets NEGBCD on BCD id)
    std::string label;             // optional; truncated to LBLEN (10) chars in output
};

// ── LineRecord ─────────────────────────────────────────────────────────────

struct LineRecord {
    std::array<int, 2 * MAXQN> qn = {}; // upper QNs then lower QNs; unused entries = 0
    int    nqn = 0;                      // number of QN pairs actually used
    double freq = 0.0;                   // observed frequency (MHz)
    double err = 1e-7;                   // measurement uncertainty (MHz)
    double weight = 1.0;                 // line weight
    std::string blend_tag;               // non-empty triggers blend-grouping in linein
};

// ── DipoleMoment ───────────────────────────────────────────────────────────

struct DipoleMoment {
    int64_t id = 0;                       // dipole identifier (decimal); builder BCD-packs it
    double  value = 0.0;                  // dipole moment (Debye)
    bool    starts_new_component = false; // true → NEGBCD flag on BCD id (new component group)
    std::string label;                    // optional label (text after '/' in .int file)
};

// ── Engine options (SPINV) ─────────────────────────────────────────────────

struct VibState {
    // Vibrational/electronic state configuration — one entry per option card.
    // vibs[0] sets defaults for all states; additional entries override per-state.
    int    index = 0;   // 0-based vibrational state index (card 1 always configures state 0)
    int    knmin = 0;   // minimum K quantum number
    int    knmax = 359; // maximum K; set both knmin=knmax=0 for a linear molecule
    // Axis for statistical weight: 1=a, 2=b, 3=c, 4=A(C2), 5=B(C2), 6=3-fold,
    // 7=A/E(C4), 8=B(C4), 9=5-fold, 10=A/E2(C6), 11=B/E1(C6). Negative → I_tot spin basis.
    int    stat_weight_axis = 1;
    int    iwtpl = 1;   // statistical weight for even (plus) parity states
    int    iwtmn = 1;   // statistical weight for odd (minus) parity states
    // vibrational symmetry (only meaningful on last VibState; builder injects -1.0 sentinel
    // for earlier states automatically — do not set to -1.0 here).
    double vsym = 0.0;
    // E-symmetry weight (EWT field). Default 99 = "ignored". Negative sign on vibs[0]
    // triggers EWTFAC=1000 mode (default EWTFAC=100). Relevant only for IAX >= 6.
    int    esym_weight = 99;
    // Spin degeneracy (2I+1) for each spin species, in units/tens/hundreds digit order.
    // Examples: [2] = one spin-1/2 nucleus; [3] = one spin-1 or one triplet electron spin;
    // [1] = spin-0 placeholder (e.g. C2v molecule with no nuclear spins).
    // For a molecule with no spin species at all, leave empty.
    std::vector<int> spin_degeneracies;
    // True if SPIND is negative: selects symmetric-top quantum numbers (K quantum number).
    // Set True for linear molecules and open-shell (electron-spin) systems.
    // False (default) for standard asymmetric-top molecules with nuclear spins only.
    bool symmetric_rotor_quanta = false;
};

struct SpinvOptions {
    // Binary flags for which inter-state couplings to include (IXX field):
    //   bit 0 set → no ΔN≠0; bit 1 → no ΔJ; bit 2 → no ΔF1; etc. Default 0 = include all.
    int    inclusion_flags = 0;
    // Eigenvalue ordering within Wang sub-blocks (DIAG field, card 1 only):
    //   0=energy, 1=full projection, 2=energy within Wang, 3=τ=Ka-Kc, 4=<Kz²>, 5=diag order.
    int    diag_order = 0;
    // Phase/operator flags (XOPT field): phase_flags = PHASE + 10*NEWLZ + 20*NOFC + 40*G12.
    //   PHASE 0-8 forces phase convention (0=none, 8=standard).
    int    phase_flags = 0;
    bool   oblate = false;    // oblate rotor: z=c, y=b, x=a (default false = prolate: z=a)
    std::string nam_file;     // parameter-label file; empty = engine default ("sping.nam")
    std::vector<VibState> vibs; // at least one entry required
};

// ── Engine options (DPI) ───────────────────────────────────────────────────

struct DpiOptions {
    int isdgn = 1;  // spin degeneracy
    int nvib  = 1;  // number of vibrational states
};

// ── Discriminated engine options ───────────────────────────────────────────

struct EngineOptions {
    enum class Kind { Spinv, Dpi };
    Kind kind = Kind::Spinv;
    SpinvOptions spinv;   // valid iff kind == Spinv
    DpiOptions   dpi;     // valid iff kind == Dpi
};

// ── CatControl ─────────────────────────────────────────────────────────────

struct CatControl {
    int    iflg = 0;
    long   itag = 999;
    double qrot = 1000.0;
    int    inblk = 0;
    int    lblk = 0;
    double thrsh = -100.0;
    double thrsh1 = -100.0;
    double fqmax = 9999.99;
    double tmq = 300.0;
    int    maxv = -1;
};

// ── Top-level fit / cat inputs ─────────────────────────────────────────────

struct FitInput {
    std::string title;
    int    n_iterations = 1;
    double marquardt_param = 0.0;
    double max_obs_calc_err = 1e6;
    double param_err_scale = 1.0;    // parfac
    double freq_scale = 1.0;         // fqfacq
    size_t max_lines = 32767;        // limlin
    bool   extended_qn = false;      // true → catqn = MAXQN instead of MAXCAT
    int    nxpar = 0;                // number of parameters excluded if index < 0
    EngineOptions engine_options;
    std::vector<Parameter>  parameters;
    std::vector<double>     variance; // packed upper triangular, length nfit*(nfit+1)/2;
                                      // empty → default (diagonals from a_priori_error)
    std::vector<LineRecord> lines;
    // Internal: raw .lin line strings from parse_fit_files; bypasses LineRecord formatting.
    // Set by parse_fit_files(); ignored if empty (build_fit_input uses lines instead).
    std::vector<std::string> raw_lines;
};

struct CatInput {
    std::string title;           // from .var/fitted.toml (Hamiltonian title)
    std::string int_title;       // from .int/dipoles.toml (catalog title for .out line 1)
    bool extended_qn = false;        // true → catqn = MAXQN instead of MAXCAT
    CatControl control;
    std::vector<DipoleMoment> dipoles;
    EngineOptions engine_options;
    std::vector<Parameter> parameters;
    std::vector<double>    variance;  // packed upper triangular (same shape as FitInput)
};

#endif // INPUT_SCHEMA_HPP
