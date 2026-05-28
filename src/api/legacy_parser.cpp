#include "api/legacy_parser.hpp"
#include "api/builders.hpp"  // for CalFitInput/CalCatInput types, LBLEN
#include "api/InputSchema.hpp"
#include "splib/ulib.h"
#include "splib/calpgm_types.h"
#include "common/CalError.hpp"

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <inttypes.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// ── Helpers ──────────────────────────────────────────────────────────────────

// Read entire file into a vector of lines (no trailing newline in each element).
static std::vector<std::string> slurp_lines(const std::string &path)
{
    FILE *f = fopen(path.c_str(), "r");
    if (!f)
        throw IoError("Cannot open file: " + path, CalErrorCode::FileOpenFailed);

    std::vector<std::string> lines;
    char buf[NDCARD];
    while (fgetstr(buf, NDCARD, f) > 0)
        lines.push_back(buf);

    fclose(f);
    return lines;
}

// Join lines into one string, each terminated by '\n'.
static std::string join_lines(const std::vector<std::string> &lines, int from, int to)
{
    std::string out;
    for (int i = from; i < to && i < (int)lines.size(); ++i) {
        out += lines[i];
        if (out.empty() || out.back() != '\n')
            out += '\n';
    }
    return out;
}

// ── SPINV option card parser ──────────────────────────────────────────────────

// Parse SPINV option cards starting at lines[offset].
// Returns SpinvOptions; sets lines_consumed to number of cards read.
static SpinvOptions parse_spinv_option_lines(const std::vector<std::string> &lines,
                                              int offset, int &lines_consumed)
{
    static const int NDECSPIN = 11;
    static const int NVEC = 11;

    SpinvOptions so;
    lines_consumed = 0;

    double rvec0[NVEC] = {1, 0, 359, 0, 1, 1, 1, 0, 99, 0, 0};
    bool first_card = true;

    for (int ci = offset; ci < (int)lines.size(); ++ci) {
        const std::string &card = lines[ci];
        const char *cstr = card.c_str();

        // Detect leading letter (ctyp) for namfil
        char ctyp = '\0';
        if (!card.empty() && isalpha((unsigned char)card[0]))
            ctyp = card[0];
        else if (card.size() > 1 && isalpha((unsigned char)card[1]))
            ctyp = card[1];

        // Parse spin integer via getbcd
        bcd_t bcdspin[11];
        bcdspin[0] = (bcd_t)NDECSPIN;
        int iend = getbcd(cstr, bcdspin, NDECSPIN);
        if (iend <= 0)
            break;  // no valid integer → end of option cards

        bool sym_rotor_quanta = (NEGBCD(bcdspin[0]) != 0);

        // Decode spin_degeneracies from bcdspin using putbcd → decimal string
        char spin_str[24] = {};
        putbcd(spin_str, sizeof(spin_str) - 1, bcdspin);

        std::vector<int> spin_degs;
        std::string spin_digits;
        for (char c : std::string(spin_str)) {
            if (c >= '0' && c <= '9') spin_digits += c;
        }
        if (spin_digits.empty() || spin_digits == "0") {
            // no spin species
        } else {
            // digits are in reversed order (highest species first): reverse back
            for (int i = (int)spin_digits.size() - 1; i >= 0; --i)
                spin_degs.push_back(spin_digits[i] - '0');
        }

        // Parse 11 floats via pcard (rvec defaults from rvec0, with vsym override to 0)
        double rvec[NVEC];
        for (int i = 0; i < NVEC; ++i) rvec[i] = rvec0[i];
        rvec[7] = 0.0;  // vsym default (not from rvec0)
        pcard(&cstr[iend], rvec, NVEC, NULL);

        int lopt   = (int)rvec[0];
        int knmin  = (int)rvec[1];
        int knmax  = (int)rvec[2];
        // rvec[3] = inclusion_flags (card 1 only)
        int stat_wt_axis = (int)rvec[4];
        int iwtpl  = (int)rvec[5];
        int iwtmn  = (int)rvec[6];
        double vsym_raw = rvec[7];  // raw value; -1.0 on non-last cards is a sentinel
        int esym_wt = (int)rvec[8];
        // rvec[9]  = diag_order (card 1 only)
        // rvec[10] = phase_flags (card 1 only)

        ++lines_consumed;

        bool is_last = (vsym_raw >= -0.5);

        VibState vib;
        vib.knmin                = knmin;
        vib.knmax                = knmax;
        vib.stat_weight_axis     = stat_wt_axis;
        vib.iwtpl                = iwtpl;
        vib.iwtmn                = iwtmn;
        // Don't propagate the -1.0 sentinel: store 0.0 for non-last states
        // (the builder injects -1.0 automatically for intermediate states).
        vib.vsym                 = is_last ? vsym_raw : 0.0;
        vib.esym_weight          = esym_wt;
        vib.spin_degeneracies    = spin_degs;
        vib.symmetric_rotor_quanta = sym_rotor_quanta;

        if (first_card) {
            first_card = false;
            int nvib = (lopt < 0) ? -lopt : lopt;
            if (nvib <= 0) nvib = 1;
            so.oblate            = (lopt < 0);
            so.inclusion_flags   = (int)rvec[3];
            so.diag_order        = (int)rvec[9];
            so.phase_flags       = (int)rvec[10];

            // Compute namfil override from ctyp (setopt: "sping.nam"[4] = ctyp)
            if (ctyp != '\0') {
                char namfil[10] = "sping.nam";
                namfil[4] = ctyp;
                so.nam_file = namfil;
            }

            vib.index = 0;
            so.vibs.resize((size_t)nvib);

            // Update rvec0 for subsequent cards (setopt does dcopy(rvec → rvec0))
            for (int i = 0; i < NVEC; ++i) rvec0[i] = rvec[i];
        } else {
            int ivib = (lopt < 0) ? -lopt : lopt;
            if (ivib >= (int)so.vibs.size())
                continue;  // out of range → skip but still consumed
            vib.index = ivib;
        }

        so.vibs[vib.index] = vib;

        // Terminator: vsym_raw >= -0.5 (last card)
        if (vsym_raw >= -0.5)
            break;
    }
    return so;
}

// ── DPI option card parser ─────────────────────────────────────────────────

static DpiOptions parse_dpi_option_line(const std::string &line)
{
    DpiOptions dpi;
    double vec[2] = {1.0, 1.0};
    pcard(line.c_str(), vec, 2, NULL);
    dpi.isdgn = (int)vec[0];
    dpi.nvib  = (int)vec[1];
    return dpi;
}

// ── Parameter line parser ─────────────────────────────────────────────────────

// Parse parameter lines starting at lines[offset].
// Returns vector<Parameter> with absolute values from the file.
// lines_consumed is set to the number of lines consumed (stops at covariance sentinel).
static std::vector<Parameter> parse_parameter_lines(const std::vector<std::string> &lines,
                                                     int offset, int &lines_consumed)
{
    std::vector<Parameter> params;
    lines_consumed = 0;
    static const double ermin = 1.0e-37;
    // erdef is the default a_priori_error used in getpar when not specified
    static const double erdef = 1.0e37;

    for (int li = offset; li < (int)lines.size(); ++li) {
        const std::string &card = lines[li];
        const char *cstr = card.c_str();

        // Parse parameter ID via getbcd (re-uses the same function as getpar)
        bcd_t idbcd[16];
        idbcd[0] = 8;  // generous nbcd (enough for any molecule type)
        int kk = getbcd(cstr, idbcd, 8);
        if (kk <= 0)
            break;  // no integer found → end of parameter section

        bool is_negbcd = (NEGBCD(idbcd[0]) != 0);

        // Decode parameter ID to int64_t via putbcd → decimal string
        char id_str[20] = {};
        putbcd(id_str, sizeof(id_str) - 1, idbcd);
        // Strip spaces, keep '-' and digits
        const char *p = id_str;
        while (*p == ' ') ++p;
        int64_t id = strtoll(p, NULL, 10);
        if (id < 0) id = -id;  // strip sign (NEGBCD flag is separate)

        // Parse value, error, weight (vec[2] is sentinel: < 1.01 → covariance)
        double vec[3] = {0.0, erdef, 2.0};
        if (pcard(&cstr[kk], vec, 3, NULL) == 0)
            break;  // empty line
        if (vec[2] < 1.01)
            break;  // reading into covariance section

        // Extract label after '/'
        std::string label;
        const char *slash = strchr(cstr, '/');
        if (slash) {
            const char *lbl = slash + 1;
            // Trim trailing whitespace
            int llen = (int)strlen(lbl);
            while (llen > 0 && (lbl[llen-1] == ' ' || lbl[llen-1] == '\n' || lbl[llen-1] == '\r'))
                --llen;
            label.assign(lbl, (size_t)llen);
            if (label.size() > (size_t)LBLEN)
                label.resize((size_t)LBLEN);
        }

        Parameter param;
        param.id    = id;
        param.value = vec[0];                              // absolute value (same as getpar's vec[0])
        param.fixed = is_negbcd;
        param.a_priori_error = is_negbcd ? -1.0 :         // NEGBCD → build_fit_input uses ermin
                               (vec[1] < ermin ? ermin : vec[1]);
        param.label = label;

        params.push_back(param);
        ++lines_consumed;
    }
    return params;
}

// ── Variance parser ──────────────────────────────────────────────────────────

// Build a minimal idpar array (ndbcd=1) for getvar: only NEGBCD flags matter.
// erpar_out gets the a_priori_error for each param (ermin for NEGBCD).
//
// Layout with ndbcd=1: param i is at idpar[i].  idpar[0] doubles as the ndbcd
// count byte (value=1) and param 0's NEGBCD byte; the legacy format requires
// param 0 to always be independent (otherwise ndbcd would be misread as 0x81),
// so NEGBCD for param 0 is always 0 in idpar[0].  Params 1..n-1 get their
// NEGBCD flags at idpar[1..n-1].
static void build_minimal_idpar(const std::vector<Parameter> &params,
                                 std::vector<unsigned char> &idpar_out,
                                 std::vector<double> &erpar_out)
{
    static const double ermin = 1.0e-37;
    size_t n = params.size();
    idpar_out.assign(n + 3, 0);
    idpar_out[0] = 1;  // ndbcd = 1; param 0 NEGBCD = 0 (always independent per format)
    erpar_out.resize(n);
    for (size_t i = 0; i < n; ++i) {
        bool is_dep = (params[i].fixed || params[i].a_priori_error < 0.0);
        // Param i's NEGBCD flag lives at idpar[i] (matches getvar's ibcd = i*ndbcd = i)
        if (i > 0)
            idpar_out[i] = is_dep ? (unsigned char)0x80 : (unsigned char)0x00;
        erpar_out[i] = is_dep ? ermin :
                       (params[i].a_priori_error < ermin ? ermin : params[i].a_priori_error);
    }
}

// Parse variance from lines[offset..end) using getvar via fmemopen.
// inpcor: 0 means try to read covariance from file; non-zero means use default.
static std::vector<double> parse_variance_lines(const std::vector<std::string> &lines,
                                                 int offset,
                                                 const std::vector<Parameter> &params,
                                                 int nfit, int inpcor)
{
    size_t n = (size_t)nfit * (nfit + 1) / 2;
    std::vector<double> var(n, 0.0);

    if (nfit <= 0) return var;

    std::vector<unsigned char> idpar;
    std::vector<double> erpar;
    build_minimal_idpar(params, idpar, erpar);

    std::string var_text = join_lines(lines, offset, (int)lines.size());
    if (var_text.empty()) {
        // No variance text: let getvar build default from erpar
        var_text = "\n";  // empty line triggers EOF in getvar
    }

    FILE *f = fmemopen((void *)var_text.data(), var_text.size(), "r");
    if (!f)
        return var;

    getvar(f, nfit, var.data(), idpar.data(), erpar.data(), inpcor);
    fclose(f);

    return var;
}

// ── parse_fit_files ──────────────────────────────────────────────────────────

FitInput parse_fit_files(const std::string &parFile,
                         const std::string &linFile,
                         EngineOptions::Kind kind)
{
    FitInput fi;

    // ── 1. Read .par file ────────────────────────────────────────────────────

    std::vector<std::string> par_lines = slurp_lines(parFile);
    if (par_lines.size() < 2)
        throw InputError("Malformed .par file (< 2 lines): " + parFile,
                         CalErrorCode::MalformedInput);

    // Title (line 0)
    {
        char buf[NDCARD];
        strncpy(buf, par_lines[0].c_str(), NDCARD - 1);
        buf[NDCARD - 1] = '\0';
        chtime(buf, 82);
        fi.title = buf;
    }

    // Header (line 1): npar, max_lines, n_iterations, nxpar, marqp0, xerrmx, parfac, fqfac
    {
        double dvec[8] = {100.0, 32767.0, 1.0, 0.0, 0.0, 1e6, 1.0, 1.0};
        pcard(par_lines[1].c_str(), dvec, 8, NULL);

        // negative max_lines signals extended QN
        fi.extended_qn = (dvec[1] < 0.0);
        if (fi.extended_qn) dvec[1] = -dvec[1];

        fi.n_iterations    = (int)dvec[2];
        fi.nxpar           = (int)dvec[3];
        fi.marquardt_param = (dvec[4] < 0.0) ? 0.0 : dvec[4];
        fi.max_obs_calc_err= dvec[5];
        fi.param_err_scale = dvec[6];
        fi.freq_scale      = dvec[7];
        fi.max_lines       = (size_t)dvec[1];
    }

    // ── 2. Parse option cards (line 2+) ─────────────────────────────────────

    int opt_consumed = 0;
    if (kind == EngineOptions::Kind::Dpi) {
        if ((int)par_lines.size() < 3)
            throw InputError("Missing DPI option card in .par file", CalErrorCode::MalformedInput);
        fi.engine_options.kind = EngineOptions::Kind::Dpi;
        fi.engine_options.dpi  = parse_dpi_option_line(par_lines[2]);
        opt_consumed = 1;
    } else {
        fi.engine_options.kind  = EngineOptions::Kind::Spinv;
        fi.engine_options.spinv = parse_spinv_option_lines(par_lines, 2, opt_consumed);
    }

    // ── 3. Parse parameter lines ─────────────────────────────────────────────

    int par_offset = 2 + opt_consumed;
    int par_consumed = 0;
    fi.parameters = parse_parameter_lines(par_lines, par_offset, par_consumed);

    // ── 4. Parse variance (if present) ───────────────────────────────────────

    // Count independent parameters
    int nfit = 0;
    for (auto &p : fi.parameters)
        if (!p.fixed && p.a_priori_error >= 0.0) ++nfit;

    // inpcor: 0 = all expected params were read → try to read covariance
    // (mirrors getpar's return value: i - n = 0 when all params read)
    int variance_offset = par_offset + par_consumed;
    int inpcor = 0;  // try to read covariance from file
    fi.variance = parse_variance_lines(par_lines, variance_offset, fi.parameters, nfit, inpcor);

    // If variance came out as default (all-zero diagonals), clear it so
    // build_fit_input will build the correct default from a_priori_error.
    // Only keep it if it has actual content (any non-zero off-diagonal).
    bool all_diag = true;
    for (int j = 0; j < nfit && all_diag; ++j) {
        for (int i = 0; i < j && all_diag; ++i) {
            size_t idx = (size_t)j*(j+1)/2 + i;
            if (fi.variance[idx] != 0.0) all_diag = false;
        }
    }
    if (all_diag)
        fi.variance.clear();  // let build_fit_input use default

    // ── 5. Read .lin file as raw strings ─────────────────────────────────────

    fi.raw_lines = slurp_lines(linFile);

    return fi;
}

// ── parse_cat_files ──────────────────────────────────────────────────────────

CatInput parse_cat_files(const std::string &varFile,
                         const std::string &intFile,
                         EngineOptions::Kind kind)
{
    CatInput ci;

    // ── 1. Read .int file ────────────────────────────────────────────────────

    std::vector<std::string> int_lines = slurp_lines(intFile);
    if (int_lines.size() < 2)
        throw InputError("Malformed .int file (< 2 lines): " + intFile,
                         CalErrorCode::MalformedInput);

    // Title (line 0)
    {
        char buf[NDCARD];
        strncpy(buf, int_lines[0].c_str(), NDCARD - 1);
        buf[NDCARD - 1] = '\0';
        chtime(buf, 82);
        ci.title = buf;
    }

    // Control line (line 1)
    {
        double dvec[10] = {0, 999, 1000, 0, 0, -100, -100, 9999.99, 300., -1};
        pcard(int_lines[1].c_str(), dvec, 10, NULL);
        ci.control.iflg   = (int)dvec[0];
        ci.control.itag   = (long)dvec[1];
        ci.control.qrot   = (dvec[2] < 1.0) ? 1.0 : dvec[2];
        ci.control.inblk  = (int)dvec[3];
        ci.control.lblk   = (int)dvec[4];
        ci.control.thrsh  = dvec[5];
        ci.control.thrsh1 = dvec[6];
        ci.control.fqmax  = dvec[7];
        ci.control.tmq    = dvec[8];
        ci.control.maxv   = (int)dvec[9];
    }

    // Dipole lines (lines 2+)
    for (int li = 2; li < (int)int_lines.size(); ++li) {
        const char *cstr = int_lines[li].c_str();
        bcd_t idbcd[8];
        idbcd[0] = NDECDIP;
        int kk = getbcd(cstr, idbcd, NDECDIP);
        if (kk <= 0) break;

        bool negbcd = (NEGBCD(idbcd[0]) != 0);
        char id_str[16] = {};
        putbcd(id_str, sizeof(id_str) - 1, idbcd);
        const char *pp = id_str;
        while (*pp == ' ') ++pp;
        int64_t id = strtoll(pp, NULL, 10);
        if (id < 0) id = -id;

        double val = 0.0;
        pcard(&cstr[kk], &val, 1, NULL);

        DipoleMoment dm;
        dm.id    = id;
        dm.value = val;
        dm.starts_new_component = negbcd && (li > 2);
        ci.dipoles.push_back(dm);
    }

    // ── 2. Read .var file ────────────────────────────────────────────────────

    std::vector<std::string> var_lines = slurp_lines(varFile);
    if (var_lines.size() < 2)
        throw InputError("Malformed .var file (< 2 lines): " + varFile,
                         CalErrorCode::MalformedInput);

    // Skip title line (var_lines[0])
    // Header (var_lines[1]): npar, catqn_flag
    {
        double dvec[2] = {0, 1};
        pcard(var_lines[1].c_str(), dvec, 2, NULL);
        ci.extended_qn = (dvec[1] < 0.0);
        // npar is informational; we read however many lines are present
    }

    // ── 3. Parse option cards (var_lines[2]+) ────────────────────────────────

    int opt_consumed = 0;
    if (kind == EngineOptions::Kind::Dpi) {
        if ((int)var_lines.size() < 3)
            throw InputError("Missing DPI option card in .var file", CalErrorCode::MalformedInput);
        ci.engine_options.kind = EngineOptions::Kind::Dpi;
        ci.engine_options.dpi  = parse_dpi_option_line(var_lines[2]);
        opt_consumed = 1;
    } else {
        ci.engine_options.kind  = EngineOptions::Kind::Spinv;
        ci.engine_options.spinv = parse_spinv_option_lines(var_lines, 2, opt_consumed);
    }

    // ── 4. Parse parameter lines ─────────────────────────────────────────────

    int par_offset = 2 + opt_consumed;
    int par_consumed = 0;
    ci.parameters = parse_parameter_lines(var_lines, par_offset, par_consumed);

    // ── 5. Parse variance from remaining .var lines ───────────────────────────

    int nfit = 0;
    for (auto &p : ci.parameters)
        if (!p.fixed && p.a_priori_error >= 0.0) ++nfit;

    int var_offset = par_offset + par_consumed;
    // inpcor=0 → try to read covariance (normal for .var files)
    ci.variance = parse_variance_lines(var_lines, var_offset, ci.parameters, nfit, 0);

    return ci;
}
