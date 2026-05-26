/**
 * toml_io.cpp — TOML file I/O for SPFIT/SPCAT
 *
 * Converts between the typed InputSchema structs and TOML files via toml++.
 * The mol.var.toml format mirrors the legacy .var file: it carries engine
 * options + fitted parameters + packed variance so spcat can read it without
 * the original .par file.
 */

#define TOML_EXCEPTIONS 0  // use error_code returns, convert to exceptions below
#include "third_party/tomlplusplus/toml.hpp"

#include "api/toml_io.hpp"
#include "api/InputSchema.hpp"
#include "spfit/CalFit.hpp"
#include "spcat/CalCat.hpp"
#include "common/CalError.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ── Internal helpers ──────────────────────────────────────────────────────

static VibState vib_from_toml(const toml::table &t)
{
    VibState v;
    if (auto x = t["knmin"].value<int64_t>())   v.knmin = (int)*x;
    if (auto x = t["knmax"].value<int64_t>())   v.knmax = (int)*x;
    if (auto x = t["stat_weight_axis"].value<int64_t>()) v.stat_weight_axis = (int)*x;
    if (auto x = t["iwtpl"].value<int64_t>())   v.iwtpl = (int)*x;
    if (auto x = t["iwtmn"].value<int64_t>())   v.iwtmn = (int)*x;
    if (auto x = t["vsym"].value<double>())     v.vsym = *x;
    if (auto x = t["esym_weight"].value<int64_t>()) v.esym_weight = (int)*x;
    if (auto x = t["symmetric_rotor_quanta"].value<bool>()) v.symmetric_rotor_quanta = *x;
    if (auto *arr = t["spin_degeneracies"].as_array()) {
        arr->for_each([&](auto &&el) {
            if constexpr (toml::is_integer<decltype(el)>)
                v.spin_degeneracies.push_back((int)*el);
        });
    }
    return v;
}

static toml::table vib_to_toml(const VibState &v)
{
    toml::table t;
    if (v.knmin != 0)   t.insert("knmin",   v.knmin);
    if (v.knmax != 359) t.insert("knmax",   v.knmax);
    if (v.stat_weight_axis != 1) t.insert("stat_weight_axis", v.stat_weight_axis);
    if (v.iwtpl != 1)   t.insert("iwtpl",   v.iwtpl);
    if (v.iwtmn != 1)   t.insert("iwtmn",   v.iwtmn);
    if (v.vsym != 0.0)  t.insert("vsym",    v.vsym);
    if (v.esym_weight != 99) t.insert("esym_weight", v.esym_weight);
    if (v.symmetric_rotor_quanta) t.insert("symmetric_rotor_quanta", true);
    if (!v.spin_degeneracies.empty()) {
        toml::array arr;
        for (int d : v.spin_degeneracies) arr.push_back(d);
        t.insert("spin_degeneracies", std::move(arr));
    }
    return t;
}

static SpinvOptions spinv_from_toml(const toml::table &t)
{
    SpinvOptions s;
    if (auto x = t["inclusion_flags"].value<int64_t>()) s.inclusion_flags = (int)*x;
    if (auto x = t["diag_order"].value<int64_t>())      s.diag_order = (int)*x;
    if (auto x = t["phase_flags"].value<int64_t>())     s.phase_flags = (int)*x;
    if (auto x = t["oblate"].value<bool>())             s.oblate = *x;
    if (auto x = t["nam_file"].value<std::string>())    s.nam_file = *x;
    if (auto *arr = t["vibs"].as_array()) {
        arr->for_each([&](auto &&el) {
            if constexpr (toml::is_table<decltype(el)>)
                s.vibs.push_back(vib_from_toml(*el.as_table()));
        });
    }
    return s;
}

static toml::table spinv_to_toml(const SpinvOptions &s)
{
    toml::table t;
    if (s.inclusion_flags != 0) t.insert("inclusion_flags", s.inclusion_flags);
    if (s.diag_order != 0)      t.insert("diag_order", s.diag_order);
    if (s.phase_flags != 0)     t.insert("phase_flags", s.phase_flags);
    if (s.oblate)               t.insert("oblate", true);
    if (!s.nam_file.empty())    t.insert("nam_file", s.nam_file);
    if (!s.vibs.empty()) {
        toml::array arr;
        for (const auto &v : s.vibs) arr.push_back(vib_to_toml(v));
        t.insert("vibs", std::move(arr));
    }
    return t;
}

static DpiOptions dpi_from_toml(const toml::table &t)
{
    DpiOptions d;
    if (auto x = t["isdgn"].value<int64_t>()) d.isdgn = (int)*x;
    if (auto x = t["nvib"].value<int64_t>())  d.nvib  = (int)*x;
    return d;
}

static toml::table dpi_to_toml(const DpiOptions &d)
{
    toml::table t;
    if (d.isdgn != 1) t.insert("isdgn", d.isdgn);
    if (d.nvib  != 1) t.insert("nvib",  d.nvib);
    return t;
}

static EngineOptions engine_options_from_toml(const toml::table &t)
{
    EngineOptions e;
    std::string kind = t["kind"].value_or(std::string("spinv"));
    if (kind == "dpi") {
        e.kind = EngineOptions::Kind::Dpi;
        if (auto *sub = t["dpi"].as_table()) e.dpi = dpi_from_toml(*sub);
    } else {
        e.kind = EngineOptions::Kind::Spinv;
        if (auto *sub = t["spinv"].as_table()) e.spinv = spinv_from_toml(*sub);
    }
    return e;
}

static toml::table engine_options_to_toml(const EngineOptions &e)
{
    toml::table t;
    if (e.kind == EngineOptions::Kind::Dpi) {
        t.insert("kind", "dpi");
        auto sub = dpi_to_toml(e.dpi);
        if (!sub.empty()) t.insert("dpi", std::move(sub));
    } else {
        t.insert("kind", "spinv");
        auto sub = spinv_to_toml(e.spinv);
        if (!sub.empty()) t.insert("spinv", std::move(sub));
    }
    return t;
}

static Parameter parameter_from_toml(const toml::table &t)
{
    Parameter p;
    if (auto x = t["id"].value<int64_t>())    p.id    = *x;
    if (auto x = t["value"].value<double>())  p.value = *x;
    if (auto x = t["a_priori_error"].value<double>()) p.a_priori_error = *x;
    if (auto x = t["fixed"].value<bool>())    p.fixed = *x;
    if (auto x = t["label"].value<std::string>()) p.label = *x;
    return p;
}

// var.toml uses 'error' (fitted sigma) in place of 'a_priori_error'
static Parameter var_parameter_from_toml(const toml::table &t)
{
    Parameter p = parameter_from_toml(t);
    if (auto x = t["error"].value<double>()) p.a_priori_error = *x;
    return p;
}

static toml::table var_parameter_to_toml(const Parameter &src, double fitted_val, double err)
{
    toml::table t;
    t.insert("id",    (int64_t)src.id);
    t.insert("value", fitted_val);
    t.insert("error", err);
    if (src.fixed) t.insert("fixed", true);
    if (!src.label.empty()) t.insert("label", src.label);
    return t;
}

static LineRecord line_record_from_toml(const toml::table &t)
{
    LineRecord lr;
    if (auto *arr = t["qn"].as_array()) {
        int i = 0;
        arr->for_each([&](auto &&el) {
            if constexpr (toml::is_integer<decltype(el)>) {
                if (i < (int)lr.qn.size()) lr.qn[(size_t)i++] = (int)*el;
            }
        });
    }
    if (auto x = t["nqn"].value<int64_t>()) lr.nqn = (int)*x;
    if (auto x = t["freq"].value<double>()) lr.freq = *x;
    if (auto x = t["err"].value<double>())  lr.err  = *x;
    if (auto x = t["weight"].value<double>()) lr.weight = *x;
    if (auto x = t["blend_tag"].value<std::string>()) lr.blend_tag = *x;
    return lr;
}

static CatControl cat_control_from_toml(const toml::table &t)
{
    CatControl c;
    if (auto x = t["iflg"].value<int64_t>())   c.iflg   = (int)*x;
    if (auto x = t["itag"].value<int64_t>())   c.itag   = (long)*x;
    if (auto x = t["qrot"].value<double>())    c.qrot   = *x;
    if (auto x = t["inblk"].value<int64_t>())  c.inblk  = (int)*x;
    if (auto x = t["lblk"].value<int64_t>())   c.lblk   = (int)*x;
    if (auto x = t["thrsh"].value<double>())   c.thrsh  = *x;
    if (auto x = t["thrsh1"].value<double>())  c.thrsh1 = *x;
    if (auto x = t["fqmax"].value<double>())   c.fqmax  = *x;
    if (auto x = t["tmq"].value<double>())     c.tmq    = *x;
    if (auto x = t["maxv"].value<int64_t>())   c.maxv   = (int)*x;
    return c;
}

static DipoleMoment dipole_from_toml(const toml::table &t)
{
    DipoleMoment dm;
    if (auto x = t["id"].value<int64_t>())    dm.id    = *x;
    if (auto x = t["value"].value<double>())  dm.value = *x;
    if (auto x = t["starts_new_component"].value<bool>()) dm.starts_new_component = *x;
    return dm;
}

static toml::parse_result parse_file(const std::string &path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        throw IoError("Cannot open TOML file: " + path, CalErrorCode::FileOpenFailed);
    std::ostringstream ss;
    ss << f.rdbuf();
    return toml::parse(ss.str());
}

static const toml::table &require_table(const toml::parse_result &r, const std::string &path)
{
    if (!r) {
        std::string msg = "TOML parse error in " + path + ": ";
        msg += std::string(r.error().description());
        throw InputError(msg, CalErrorCode::MalformedInput);
    }
    return r.table();
}

// ── Public API ────────────────────────────────────────────────────────────

FitInput load_fit_input_toml(const std::string &path)
{
    auto result = parse_file(path);
    const toml::table &tbl = require_table(result, path);

    FitInput fi;
    fi.title            = tbl["title"].value_or(std::string{});
    fi.n_iterations     = (int)tbl["n_iterations"].value_or((int64_t)1);
    fi.marquardt_param  = tbl["marquardt_param"].value_or(0.0);
    fi.max_obs_calc_err = tbl["max_obs_calc_err"].value_or(1e6);
    fi.param_err_scale  = tbl["param_err_scale"].value_or(1.0);
    fi.freq_scale       = tbl["freq_scale"].value_or(1.0);
    fi.max_lines        = (size_t)tbl["max_lines"].value_or((int64_t)32767);
    fi.extended_qn      = tbl["extended_qn"].value_or(false);
    fi.nxpar            = (int)tbl["nxpar"].value_or((int64_t)0);

    if (auto *et = tbl["engine_options"].as_table())
        fi.engine_options = engine_options_from_toml(*et);

    if (auto *pa = tbl["parameters"].as_array()) {
        pa->for_each([&](auto &&el) {
            if constexpr (toml::is_table<decltype(el)>)
                fi.parameters.push_back(parameter_from_toml(*el.as_table()));
        });
    }
    if (auto *va = tbl["variance"].as_array()) {
        va->for_each([&](auto &&el) {
            if constexpr (toml::is_floating_point<decltype(el)>)
                fi.variance.push_back(*el);
        });
    }
    if (auto *la = tbl["lines"].as_array()) {
        la->for_each([&](auto &&el) {
            if constexpr (toml::is_table<decltype(el)>)
                fi.lines.push_back(line_record_from_toml(*el.as_table()));
        });
    }
    return fi;
}

CatInput load_cat_input_toml(const std::string &var_path, const std::string &int_path)
{
    auto var_result = parse_file(var_path);
    const toml::table &var_tbl = require_table(var_result, var_path);
    auto int_result = parse_file(int_path);
    const toml::table &int_tbl = require_table(int_result, int_path);

    CatInput ci;
    ci.title = var_tbl["title"].value_or(std::string{});

    if (auto *et = var_tbl["engine_options"].as_table())
        ci.engine_options = engine_options_from_toml(*et);

    if (auto *pa = var_tbl["parameters"].as_array()) {
        pa->for_each([&](auto &&el) {
            if constexpr (toml::is_table<decltype(el)>)
                ci.parameters.push_back(var_parameter_from_toml(*el.as_table()));
        });
    }
    if (auto *va = var_tbl["variance"].as_array()) {
        va->for_each([&](auto &&el) {
            if constexpr (toml::is_floating_point<decltype(el)>)
                ci.variance.push_back(*el);
        });
    }

    if (auto *ct = int_tbl["control"].as_table())
        ci.control = cat_control_from_toml(*ct);
    if (auto *da = int_tbl["dipoles"].as_array()) {
        da->for_each([&](auto &&el) {
            if constexpr (toml::is_table<decltype(el)>)
                ci.dipoles.push_back(dipole_from_toml(*el.as_table()));
        });
    }
    return ci;
}

void save_fit_output_toml(const CalFitOutput &out, const FitInput &src_fi,
                          const std::string &path)
{
    toml::table tbl;
    tbl.insert("title",   src_fi.title);
    tbl.insert("xsqbest", out.xsqbest);
    tbl.insert("itr",     (int64_t)out.itr);
    tbl.insert("engine_options", engine_options_to_toml(src_fi.engine_options));

    toml::array params_arr;
    size_t n = std::min({src_fi.parameters.size(), out.par.size(), out.erpar.size()});
    for (size_t i = 0; i < n; ++i)
        params_arr.push_back(var_parameter_to_toml(src_fi.parameters[i],
                                                   out.par[i], out.erpar[i]));
    tbl.insert("parameters", std::move(params_arr));

    if (!out.var_final_for_output.empty()) {
        toml::array var_arr;
        for (double v : out.var_final_for_output) var_arr.push_back(v);
        tbl.insert("variance", std::move(var_arr));
    }

    std::ofstream f(path, std::ios::binary);
    if (!f)
        throw IoError("Cannot write TOML file: " + path, CalErrorCode::FileOpenFailed);
    f << tbl;
}

void save_cat_output_toml(const CalCatOutput &out, const std::string &path)
{
    toml::table tbl;
    tbl.insert("nline", (int64_t)out.nline);

    toml::array temp_arr;
    for (int i = 0; i < out.ntemp; ++i) temp_arr.push_back(out.temp[i]);
    tbl.insert("temp", std::move(temp_arr));

    toml::array qsum_arr;
    for (int i = 0; i < out.ntemp; ++i) qsum_arr.push_back(out.qsum[i]);
    tbl.insert("qsum", std::move(qsum_arr));

    toml::array cat_arr;
    for (const auto &line : out.cat_lines) cat_arr.push_back(line);
    tbl.insert("cat_lines", std::move(cat_arr));

    std::ofstream f(path, std::ios::binary);
    if (!f)
        throw IoError("Cannot write TOML file: " + path, CalErrorCode::FileOpenFailed);
    f << tbl;
}
