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
};

// ── Engine options (SPINV) ─────────────────────────────────────────────────

struct VibState {
    // Vibrational state configuration — one entry per option card after the first.
    // The first entry (vibs[0]) also carries global setup for card 1.
    int    index = 0;   // 0-based vibrational state index; card 1 always configures index 0
    int    knmin = 0;   // minimum K quantum number
    int    knmax = 359; // maximum K (clamped to MAXN_DIRCOS = 359); 0 → linear molecule
    int    iax = 1;     // axis / symmetry selector (0–11; negative enables special knmin treatment)
    int    iwtpl = 1;   // statistical weight (plus parity)
    int    iwtmn = 1;   // statistical weight (minus parity)
    double vsym = 0.0;  // vibrational symmetry parameter (builder auto-sets intermediate cards)
    int    ewt0 = 99;   // energy/weight flag; negative value on vibs[0] sets esymdec = 1000
    // Nuclear spin encoding: 2*I values (e.g. 0 = none, 1 = spin-1/2, 2 = spin-1)
    // Maximum 9 spins; each 2*I must be in [1, 9] (isbig=0 BCD format).
    // For "isbig" spins (2*I > 9) a ValidationError is thrown.
    std::vector<int> nuclear_spins;
    // NEGBCD flag on the spin BCD field: required for linear molecules (sets nqnn=1).
    // Parser sets this automatically; users of linear molecules must set it explicitly.
    bool negbcd_spin = false;
};

struct SpinvOptions {
    int    ixz = 0;           // x-z plane selection (card 1 field rvec[3])
    int    idiag = 0;         // diagonalisation option (card 1 field rvec[9])
    int    phase_flags = 0;   // packed phase: stdphase = phase_flags%10,
                              //   newlz = (phase_flags/10)%2, nofc = (phase_flags/10/2)%2,
                              //   g12 = (phase_flags/10/4)%2 (see spinv_setup.cpp:620–632)
    bool   oblate = false;    // global oblate flag: negates lopt on card 1
    std::string nam_file;     // parameter-name file path; empty = engine default ("sping.nam")
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
    std::string title;
    CatControl control;
    std::vector<DipoleMoment> dipoles;
    EngineOptions engine_options;
    std::vector<Parameter> parameters;
    std::vector<double>    variance;  // packed upper triangular (same shape as FitInput)
};

#endif // INPUT_SCHEMA_HPP
