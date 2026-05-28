#include "api/builders.hpp"
#include "api/InputSchema.hpp"
#include "splib/ulib.h"
#include "splib/calpgm_types.h"
#include "splib/lsqfit.h"
#include "common/CalError.hpp"

#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <inttypes.h>
#include <sstream>

// ── Internal helpers ─────────────────────────────────────────────────────────

// Encode a decimal parameter/dipole ID into BCD, with optional NEGBCD.
static void encode_id_to_bcd(int64_t id, bool negbcd, int ndbcd, bcd_t *out)
{
    char id_str[64];
    if (negbcd)
        snprintf(id_str, sizeof(id_str), "-%" PRId64, (int64_t)llabs(id));
    else
        snprintf(id_str, sizeof(id_str), "%" PRId64, id);
    out[0] = (bcd_t)ndbcd; // set count byte first (getbcd overwrites it)
    getbcd(id_str, out, ndbcd);
}

// Format a LineRecord as a fixed-width card string that getlin_str can parse.
// nqn is the number of QN pairs (upper + lower), from eng.setfmt().
static std::string format_line_record(const LineRecord &lr, int nqn)
{
    // Column layout: 2*nqn QN values in 3-char fields, then freq err weight.
    int nn = (nqn <= 6) ? 12 : 2 * nqn;
    int nc = 3 * nn; // QN section width

    std::ostringstream oss;
    // QN fields: right-justified in 3-char columns
    for (int i = 0; i < nn; ++i) {
        int qn_val = 0;
        if (i < 2 * MAXQN && i < 2 * lr.nqn)
            qn_val = lr.qn[i];
        oss << std::setw(3) << qn_val;
    }
    (void)nc; // nc == oss.str().size() at this point

    // Freq, err, weight as free-format floats (pcard is format-agnostic)
    oss << std::setprecision(10) << std::fixed << std::setw(14) << lr.freq
        << "  " << std::scientific << std::setprecision(4) << std::setw(12) << lr.err
        << "  " << std::fixed << std::setprecision(4) << lr.weight;

    return oss.str();
}

// Build default variance: packed upper-triangular with diagonals = a_priori_error
// of each independent parameter (matching getvar's "n==0" default branch).
static std::vector<double> build_default_variance(const std::vector<double> &erp_initial,
                                                   int nfit,
                                                   const std::vector<unsigned char> &idpar_data,
                                                   int ndbcd)
{
    std::vector<double> var(((size_t)nfit * (nfit + 1)) / 2, 0.0);
    int j_fit = 0;
    int ibcd = 0; // param i is at offset i*ndbcd (matches getpar/CalFit convention)
    for (int i = 0; i < (int)erp_initial.size() && j_fit < nfit; ++i, ibcd += ndbcd) {
        if (NEGBCD(idpar_data[ibcd]) == 0) { // independent
            size_t diag = (size_t)j_fit * (j_fit + 1) / 2 + j_fit;
            var[diag] = erp_initial[i];
            ++j_fit;
        }
    }
    return var;
}

// ── build_fit_input ──────────────────────────────────────────────────────────

CalFitInput build_fit_input(const FitInput &fi, CalculationEngine &eng, Logger &/*logger*/)
{
    // ── 1. Validate basic schema constraints ─────────────────────────────────

    if (fi.parameters.empty())
        throw InputError("FitInput.parameters is empty", CalErrorCode::InvalidParameter);
    if (fi.lines.empty() && fi.raw_lines.empty())
        throw InputError("FitInput.lines and FitInput.raw_lines are both empty",
                         CalErrorCode::InvalidParameter);

    int npar = (int)fi.parameters.size();
    for (int i = 0; i < npar; ++i) {
        const Parameter &p = fi.parameters[i];
        if ((int)p.label.size() > LBLEN)
            throw ValidationError("Parameter label too long (max " + std::to_string(LBLEN) + " chars)",
                                  CalErrorCode::InvalidParameter);
    }
    if (!fi.variance.empty()) {
        // count independent params to determine expected variance size
        int nfit_expected = 0;
        for (auto &p : fi.parameters)
            if (!p.fixed && !(p.a_priori_error < 0.0)) ++nfit_expected;
        size_t expected = (size_t)nfit_expected * (nfit_expected + 1) / 2;
        if (fi.variance.size() != expected)
            throw ValidationError("FitInput.variance size mismatch: expected " +
                                  std::to_string(expected) + ", got " +
                                  std::to_string(fi.variance.size()),
                                  CalErrorCode::SizeMismatch);
    }
    for (auto &lr : fi.lines) {
        if (lr.nqn > MAXQN)
            throw ValidationError("LineRecord.nqn exceeds MAXQN=" + std::to_string(MAXQN),
                                  CalErrorCode::InvalidParameter);
        if (!std::isfinite(lr.freq))
            throw ValidationError("LineRecord.freq is not finite", CalErrorCode::InvalidParameter);
    }

    // ── 2. Apply engine options ───────────────────────────────────────────────

    CalFitInput input;
    input.title   = fi.title;
    input.npar    = npar;
    input.limlin  = fi.max_lines;
    input.nitr    = fi.n_iterations;
    input.nxpar_from_file = fi.nxpar;
    input.marqp0  = (fi.marquardt_param < 0.0) ? 0.0 : fi.marquardt_param;
    input.xerrmx  = fi.max_obs_calc_err;
    input.parfac_initial = fi.param_err_scale;
    input.fqfacq  = fi.freq_scale;
    input.catqn   = fi.extended_qn ? MAXQN : MAXCAT;

    // Sanity check on limlin / npar product (mirrors CalFitIO check)
    {
        size_t max_bytes = 32ULL << 30;
        if (const char *env = getenv("PICKETT_MAX_FIT_BYTES")) {
            char *end = nullptr;
            unsigned long long v = strtoull(env, &end, 10);
            if (end != env && *end == '\0' && v > 0) max_bytes = (size_t)v;
        }
        size_t per_line = sizeof(SXLINE) + (size_t)npar * sizeof(double);
        if (per_line > 0 && input.limlin > max_bytes / per_line)
            throw ValidationError("FitInput: limlin/npar would exceed memory sanity limit",
                                  CalErrorCode::InvalidParameter);
    }

    int temp_nfmt = input.catqn;
    int temp_itd  = 2;
    int temp_ndbcd = 1;
    std::string namfil;

    input.noptn_read_by_setopt = eng.apply_options(fi.engine_options,
                                                   &temp_nfmt, &temp_itd, &temp_ndbcd,
                                                   namfil);
    if (input.noptn_read_by_setopt < 0)
        throw InputError("apply_options failed (engine returned negative count)",
                         CalErrorCode::MalformedInput);

    input.nfmt_cat_from_setopt = temp_nfmt;
    input.itd_from_setopt      = temp_itd;
    input.ndbcd_from_setopt    = temp_ndbcd;
    input.namfil_from_setopt   = namfil;

    // ── 3. Build idpar_data, par_initial, erp_initial, parlbl_data_flat ──────

    int ndbcd = temp_ndbcd;
    // idpar_data layout: param i at offset i*ndbcd (same as getpar/CalFit convention);
    // idpar[0] is both count byte and param 0's sign byte. Extra ndbcd+3 sentinel bytes.
    size_t idpar_count = (size_t)npar * ndbcd + ndbcd + 3;
    input.idpar_data.assign(idpar_count, 0);

    input.par_initial.assign(npar, 0.0);
    input.erp_initial.assign(npar, 0.0);
    input.parlbl_data_flat.assign((size_t)LBLEN * npar + 1, '\0');

    const double ermin = 1.0e-37;
    double parbase = 1.0;
    int nfit = 0;
    input.inpcor = 0;

    for (int i = 0; i < npar; ++i) {
        const Parameter &p = fi.parameters[i];
        bool is_dep = (p.fixed || p.a_priori_error < 0.0);

        // Encode BCD id — layout matches getpar: parameter i at offset i*ndbcd,
        // with idpar[0] being both the count byte and param 0's sign byte.
        bcd_t *id_slot = input.idpar_data.data() + (size_t)i * ndbcd;
        encode_id_to_bcd(p.id, is_dep, ndbcd, id_slot);

        // Parameter value
        if (!is_dep) {
            input.par_initial[i] = p.value;
            input.erp_initial[i] = (p.a_priori_error < ermin) ? ermin : p.a_priori_error;
            parbase = p.value;
            if (std::fabs(parbase) < ermin) parbase = ermin;
            ++nfit;
        } else {
            // Dependent: store ratio relative to preceding independent value
            input.par_initial[i] = p.value / parbase;
            input.erp_initial[i] = ermin;
        }

        // Label (LBLEN chars, null-padded)
        char *lbl_slot = input.parlbl_data_flat.data() + (size_t)i * LBLEN;
        if (!p.label.empty()) {
            size_t len = std::min(p.label.size(), (size_t)LBLEN);
            memcpy(lbl_slot, p.label.c_str(), len);
        }
    }
    input.nfit = nfit;

    // ── 4. Build variance matrix ──────────────────────────────────────────────

    if (!fi.variance.empty()) {
        input.var_initial_from_getvar = fi.variance;
        input.inpcor = 1;
    } else {
        input.var_initial_from_getvar = build_default_variance(
            input.erp_initial, nfit, input.idpar_data, ndbcd);
        input.inpcor = 0;
    }

    // ── 5. Format line records as legacy card strings ─────────────────────────

    if (!fi.raw_lines.empty()) {
        // Fast path from parse_fit_files: raw strings are already in the correct format
        // (setfmt not called here — avoid disturbing engine static state before CalFit::run)
        input.lineData_raw = fi.raw_lines;
    } else {
        // Get quantum number format for line record formatting
        int iqnfmt_for_lines = 0;
        int nqn = eng.setfmt(&iqnfmt_for_lines, 1);
        if (nqn <= 0) nqn = 1;
        input.lineData_raw.reserve(fi.lines.size());
        for (auto &lr : fi.lines)
            input.lineData_raw.push_back(format_line_record(lr, nqn));
    }

    return input;
}

// ── build_cat_input ──────────────────────────────────────────────────────────

CalCatInput build_cat_input(const CatInput &ci, CalculationEngine &eng, Logger &/*logger*/)
{
    // ── 1. Validate ───────────────────────────────────────────────────────────

    if (ci.parameters.empty())
        throw InputError("CatInput.parameters is empty", CalErrorCode::InvalidParameter);

    int npar = (int)ci.parameters.size();

    if (!ci.variance.empty()) {
        int nfit_expected = 0;
        for (auto &p : ci.parameters)
            if (!p.fixed && !(p.a_priori_error < 0.0)) ++nfit_expected;
        size_t expected = (size_t)nfit_expected * (nfit_expected + 1) / 2;
        if (ci.variance.size() != expected)
            throw ValidationError("CatInput.variance size mismatch: expected " +
                                  std::to_string(expected) + ", got " +
                                  std::to_string(ci.variance.size()),
                                  CalErrorCode::SizeMismatch);
    }

    // ── 2. Apply engine options ───────────────────────────────────────────────

    CalCatInput input;
    const CatControl &cc = ci.control;

    input.title  = ci.title;
    input.iflg   = cc.iflg;
    input.itag   = cc.itag;
    input.qrot   = (cc.qrot < 1.0) ? 1.0 : cc.qrot;
    input.inblk  = cc.inblk;
    input.lblk   = cc.lblk;
    input.thrsh  = cc.thrsh;
    input.thrsh1 = cc.thrsh1;
    input.fqmax  = cc.fqmax;
    input.tmq    = cc.tmq;
    input.maxv   = cc.maxv;
    input.npar   = npar;
    input.catqn  = ci.extended_qn ? MAXQN : MAXCAT;

    int temp_nfmt  = input.catqn;
    int temp_itd   = 2;
    int temp_ndbcd = 1;
    std::string namfil;

    int n_opts = eng.apply_options(ci.engine_options, &temp_nfmt, &temp_itd, &temp_ndbcd, namfil);
    if (n_opts < 0)
        throw InputError("apply_options failed for CatInput", CalErrorCode::MalformedInput);

    input.nfmt  = temp_nfmt;
    input.itd   = temp_itd;
    input.ndbcd = temp_ndbcd;

    // setfmt to populate iqnfmtv
    input.iqnfmtv.assign((size_t)temp_nfmt * 2, 0);
    input.nqn = eng.setfmt(input.iqnfmtv.data(), temp_nfmt);

    // ── 3. Build dipole data ──────────────────────────────────────────────────

    int ndip = (int)ci.dipoles.size();
    input.ndip = ndip;
    input.dip.resize(ndip, 0.0);
    input.nvdip.resize(ndip, 0);
    input.isimag.resize(ndip, -1); // -1 = unset; filled later by setint
    input.idip.resize((size_t)ndip * NDECDIP, 0);
    input.idip[0] = (bcd_t)NDECDIP;

    int k = -1;
    for (int j = 0; j < ndip; ++j) {
        const DipoleMoment &dm = ci.dipoles[j];
        bcd_t *idip_slot = input.idip.data() + (size_t)j * NDECDIP;
        encode_id_to_bcd(dm.id, dm.starts_new_component && j > 0, NDECDIP, idip_slot);
        input.dip[j] = dm.value;
        input.isimag[j] = -1;

        if (NEGBCD(input.idip[j * NDECDIP]) == 0 || j == 0)
            input.nvdip[++k] = 1;
        else
            input.nvdip[k] += 1;
    }

    // Compute npdip from iflg (replicating CalCatIO logic)
    input.npdip = 1;
    if (input.iflg >= 10 && (input.iflg % 100) >= 10) {
        if ((input.iflg % 100) >= 20)
            input.npdip = -(k + 1);
        else
            input.npdip = 1;
    }
    if (input.npdip > 0) {
        if (ndip > 0) input.nvdip[0] = ndip;
        input.npdip = 1;
    } else {
        input.npdip = k + 1;
    }

    // ── 4. Build parameter arrays ─────────────────────────────────────────────

    int ndbcd = temp_ndbcd;
    size_t idpar_count = (size_t)npar * ndbcd + ndbcd + 3;
    input.idpar.assign(idpar_count, 0);

    input.par.assign(npar, 0.0);
    input.derv.assign(npar, 0.0); // erpar for CalCat

    const double ermin = 1.0e-37;
    double parbase = 1.0;
    int nfit = 0;

    for (int i = 0; i < npar; ++i) {
        const Parameter &p = ci.parameters[i];
        bool is_dep = (p.fixed || p.a_priori_error < 0.0);

        // param i at offset i*ndbcd (same as getvar/CalCat convention)
        bcd_t *id_slot = input.idpar.data() + (size_t)i * ndbcd;
        encode_id_to_bcd(p.id, is_dep, ndbcd, id_slot);

        if (!is_dep) {
            input.par[i]  = p.value;
            input.derv[i] = (p.a_priori_error < ermin) ? ermin : p.a_priori_error;
            parbase = p.value;
            if (std::fabs(parbase) < ermin) parbase = ermin;
            ++nfit;
        } else {
            input.par[i]  = p.value / parbase;
            input.derv[i] = ermin;
        }
    }
    input.nfit = nfit;

    // ── 5. Build variance matrix ──────────────────────────────────────────────

    if (!ci.variance.empty()) {
        input.var = ci.variance;
    } else {
        // Default: diagonals from derv (a_priori_error), off-diagonals zero
        std::vector<double> var(((size_t)nfit * (nfit + 1)) / 2, 0.0);
        int j_fit = 0;
        for (int i = 0; i < npar && j_fit < nfit; ++i) {
            bcd_t *id_slot = input.idpar.data() + (size_t)i * ndbcd;
            if (NEGBCD(id_slot[0]) == 0) {
                size_t diag = (size_t)j_fit * (j_fit + 1) / 2 + j_fit;
                var[diag] = input.derv[i];
                ++j_fit;
            }
        }
        input.var = std::move(var);
    }

    return input;
}
