/**
 * bindings.cpp — nanobind module exposing SPFIT/SPCAT as the Python package _pickett.
 *
 * Session objects (FitSession / CatSession) bundle the engine lifecycle: the same
 * engine instance configured by CalFitIO/CalCatIO::readInput is passed to CalFit/CalCat,
 * which move-owns it.  Each session can therefore only be run once.
 *
 * High-level convenience wrappers (fit_from_files / cat_from_files) are the easiest
 * entry point and are re-exported by pickett/__init__.py.
 */

#include <nanobind/nanobind.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <sstream>
#include <iomanip>

#include "common/CalError.hpp"
#include "common/Logger.hpp"
#include "engine/SpinvEngine.hpp"
#include "engine/DpiEngine.hpp"
#include "spfit/CalFit.hpp"
#include "spfit/CalFitIO.hpp"
#include "spcat/CalCat.hpp"
#include "spcat/CalCatIO.hpp"
#include "spcat/OutputSink.hpp"
#include "api/InputSchema.hpp"
#include "api/builders.hpp"
#include "api/legacy_parser.hpp"

namespace nb = nanobind;
using namespace nb::literals;

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

static std::unique_ptr<CalculationEngine> make_engine(const std::string &name)
{
    if (name == "dpi")
        return std::make_unique<DpiEngine>();
    if (name != "spinv")
        throw ValidationError("Unknown engine '" + name + "'; use 'spinv' or 'dpi'",
                              CalErrorCode::InvalidParameter);
    return std::make_unique<SpinvEngine>();
}

// ---------------------------------------------------------------------------
// Session wrappers — pair a configured engine with its read input so that
// both can be forwarded into CalFit / CalCat in a single run() call.
// ---------------------------------------------------------------------------

struct FitSession {
    CalFitInput input;
    std::unique_ptr<CalculationEngine> engine;
    bool used = false;

    FitSession(const std::string &par_file,
               const std::string &lin_file,
               const std::string &engine_name)
    {
        engine = make_engine(engine_name);
        if (!CalFitIO::readInput(par_file, lin_file, input, engine, stdout))
            throw IoError("Failed to read fit input: " + par_file);
    }

    explicit FitSession(const FitInput &fi) {
        engine = make_engine(fi.engine_options.kind == EngineOptions::Kind::Dpi ? "dpi" : "spinv");
        input = build_fit_input(fi, *engine);
    }

    CalFitOutput run()
    {
        if (used)
            throw ValidationError("FitSession.run() may only be called once",
                                  CalErrorCode::InvalidParameter);
        used = true;
        CalFit fit(engine, stdout);
        CalFitOutput output;
        fit.run(input, output);
        return output;
    }
};

struct CatSession {
    CalCatInput input;
    std::unique_ptr<CalculationEngine> engine;
    FileSink out_sink; // diagnostic → stdout
    bool used = false;

    CatSession(const std::string &int_file,
               const std::string &var_file,
               const std::string &engine_name)
        : out_sink(stdout)
    {
        engine = make_engine(engine_name);
        if (!CalCatIO::readInput(int_file, var_file, input, engine, &out_sink))
            throw IoError("Failed to read cat input: " + int_file);
    }

    explicit CatSession(const CatInput &ci)
        : out_sink(stdout)
    {
        engine = make_engine(ci.engine_options.kind == EngineOptions::Kind::Dpi ? "dpi" : "spinv");
        input = build_cat_input(ci, *engine);
    }

    CalCatOutput run()
    {
        if (used)
            throw ValidationError("CatSession.run() may only be called once",
                                  CalErrorCode::InvalidParameter);
        used = true;

        MemorySink cat_sink, egy_sink, str_sink;
        CalCat cat(engine, &out_sink, &cat_sink, &egy_sink, &str_sink);

        CalCatOutput output;
        cat.run(input, output);

        output.cat_lines = cat_sink.drain_lines();
        output.egy_lines = egy_sink.drain_lines();
        output.str_lines = str_sink.drain_lines();
        output.sort_cat_lines();
        return output;
    }
};

// ---------------------------------------------------------------------------
// Convenience free functions (used by pickett/__init__.py)
// ---------------------------------------------------------------------------

static CalFitOutput fit_from_files(const std::string &par_file,
                                   const std::string &lin_file,
                                   const std::string &engine_name)
{
    FitSession s(par_file, lin_file, engine_name);
    return s.run();
}

static CalCatOutput cat_from_files(const std::string &int_file,
                                   const std::string &var_file,
                                   const std::string &engine_name)
{
    CatSession s(int_file, var_file, engine_name);
    return s.run();
}

// ---------------------------------------------------------------------------
// Module definition
// ---------------------------------------------------------------------------

NB_MODULE(_pickett, m)
{
    m.doc() = "SPFIT/SPCAT spectroscopy software — Python bindings";

    // ---- Exception hierarchy ----
    auto exc_cal     = nb::exception<CalError>      (m, "CalError");
    nb::exception<IoError>         (m, "IoError",         exc_cal.ptr());
    nb::exception<InputError>      (m, "InputError",      exc_cal.ptr());
    nb::exception<ValidationError> (m, "ValidationError", exc_cal.ptr());
    nb::exception<NumericError>    (m, "NumericError",     exc_cal.ptr());

    // ---- CalFitOutput ----
    nb::class_<CalFitOutput>(m, "CalFitOutput")
        .def_ro("par",      &CalFitOutput::par,     "Fitted parameter values")
        .def_ro("erpar",    &CalFitOutput::erpar,   "Estimated parameter errors")
        .def_ro("xsqbest",  &CalFitOutput::xsqbest, "Best RMS (obs-calc)/err after fitting")
        .def_ro("itr",      &CalFitOutput::itr,     "Number of iterations performed")
        .def_ro("variance", &CalFitOutput::var_final_for_output,
                "Packed upper-triangular variance matrix (nfit*(nfit+1)/2 elements); "
                "pass to save_fit_output so spcat gets accurate ERR values")
        .def("__repr__", [](const CalFitOutput &o) {
            std::ostringstream oss;
            oss << "CalFitOutput(xsqbest=" << o.xsqbest
                << ", itr=" << o.itr
                << ", parameters=" << o.par.size() << ")";
            return oss.str();
        });

    // ---- CalCatOutput ----
    nb::class_<CalCatOutput>(m, "CalCatOutput")
        .def_ro("nline",     &CalCatOutput::nline,     "Total catalog lines generated")
        .def_ro("cat_lines", &CalCatOutput::cat_lines, "Catalog lines (sorted by frequency)")
        .def_ro("egy_lines", &CalCatOutput::egy_lines, "Energy level lines")
        .def_ro("str_lines", &CalCatOutput::str_lines, "Intensity/strength lines")
        .def_ro("ntemp",     &CalCatOutput::ntemp,     "Number of temperature points")
        .def_prop_ro("temp", [](const CalCatOutput &o) {
            return std::vector<double>(o.temp, o.temp + o.ntemp);
        }, "Temperature values (K)")
        .def_prop_ro("qsum", [](const CalCatOutput &o) {
            return std::vector<double>(o.qsum, o.qsum + o.ntemp);
        }, "Partition function Q at each temperature")
        .def("__repr__", [](const CalCatOutput &o) {
            std::ostringstream oss;
            oss << "CalCatOutput(nline=" << o.nline
                << ", cat_lines=" << o.cat_lines.size()
                << ", ntemp=" << o.ntemp << ")";
            return oss.str();
        });

    // ---- Low-level session objects ----
    nb::class_<FitSession>(m, "FitSession",
        "SPFIT session — configure via files or FitInput struct.  Each instance may call run() once.")
        .def(nb::init<const std::string &, const std::string &, const std::string &>(),
             "par_file"_a, "lin_file"_a, "engine"_a = "spinv",
             "Open par_file and lin_file; configure the engine via setopt.")
        .def_static("from_input", [](const FitInput &fi) { return FitSession(fi); },
             "fit_input"_a, "Build a FitSession from a FitInput struct (no file I/O).")
        .def("run", &FitSession::run,
             "Run the fitting and return a CalFitOutput.  May only be called once.");

    nb::class_<CatSession>(m, "CatSession",
        "SPCAT session — configure via files or CatInput struct.  Each instance may call run() once.")
        .def(nb::init<const std::string &, const std::string &, const std::string &>(),
             "int_file"_a, "var_file"_a, "engine"_a = "spinv",
             "Open int_file and var_file; configure the engine via setopt.")
        .def_static("from_input", [](const CatInput &ci) { return CatSession(ci); },
             "cat_input"_a, "Build a CatSession from a CatInput struct (no file I/O).")
        .def("run", &CatSession::run,
             "Generate catalog and return a CalCatOutput.  May only be called once.");

    // ---- High-level convenience functions ----
    m.def("fit_from_files", &fit_from_files,
          "par_file"_a, "lin_file"_a, "engine"_a = "spinv",
          "Read SPFIT input files, run the fitter, and return CalFitOutput.");

    m.def("cat_from_files", &cat_from_files,
          "int_file"_a, "var_file"_a, "engine"_a = "spinv",
          "Read SPCAT input files, generate the catalog, and return CalCatOutput.");

    // ---- InputSchema — typed input structs ----
    //
    // Every input struct supports:
    //   (a) default construction:          x = Klass()
    //   (b) keyword-argument construction: x = Klass(field1=..., field2=...)
    //   (c) attribute assignment:          x.field = ...
    //   (d) a useful __repr__ for debugging
    //
    // All __init__ lambdas use placement-new of a default-constructed instance
    // followed by field assignment.  This lets nanobind generate per-field
    // default arguments without requiring corresponding C++ constructors.

    nb::enum_<EngineOptions::Kind>(m, "EngineKind")
        .value("Spinv", EngineOptions::Kind::Spinv)
        .value("Dpi",   EngineOptions::Kind::Dpi);

    // Small repr helpers
    // Precision 12 preserves typical spectroscopic values exactly (57635.968040505
    // round-trips) without the trailing-noise digits of the full 17-digit form.
    auto repr_double = [](double v) {
        std::ostringstream oss;
        oss << std::setprecision(12) << v;
        return oss.str();
    };
    auto repr_string = [](const std::string &s) {
        std::ostringstream oss;
        oss << '\'' << s << '\'';
        return oss.str();
    };

    nb::class_<Parameter>(m, "Parameter")
        .def("__init__", [](Parameter *self, int64_t id, double value,
                            double a_priori_error, bool fixed, std::string label) {
                 new (self) Parameter();
                 self->id = id;
                 self->value = value;
                 self->a_priori_error = a_priori_error;
                 self->fixed = fixed;
                 self->label = std::move(label);
             },
             "id"_a = 0, "value"_a = 0.0, "a_priori_error"_a = 1e37,
             "fixed"_a = false, "label"_a = "")
        .def_rw("id",             &Parameter::id,             "Decimal parameter ID")
        .def_rw("value",          &Parameter::value,          "Absolute value (MHz or cm-1)")
        .def_rw("a_priori_error", &Parameter::a_priori_error, "1-sigma uncertainty; <=1e-37 fixes parameter")
        .def_rw("fixed",          &Parameter::fixed,          "Excluded from fit if true")
        .def_rw("label",          &Parameter::label,          "Optional label (truncated to 10 chars)")
        .def("__repr__", [repr_double, repr_string](const Parameter &p) {
            std::ostringstream oss;
            oss << "Parameter(id=" << p.id
                << ", value=" << repr_double(p.value)
                << ", a_priori_error=" << repr_double(p.a_priori_error)
                << ", fixed=" << (p.fixed ? "True" : "False")
                << ", label=" << repr_string(p.label) << ")";
            return oss.str();
        });

    nb::class_<LineRecord>(m, "LineRecord")
        .def("__init__", [](LineRecord *self, const std::vector<int> &qn, int nqn,
                            double freq, double err, double weight, std::string blend_tag) {
                 new (self) LineRecord();
                 if (qn.size() > self->qn.size())
                     throw nb::value_error("qn list too long (max 2*MAXQN entries)");
                 std::copy(qn.begin(), qn.end(), self->qn.begin());
                 self->nqn = nqn;
                 self->freq = freq;
                 self->err = err;
                 self->weight = weight;
                 self->blend_tag = std::move(blend_tag);
             },
             "qn"_a = std::vector<int>{}, "nqn"_a = 0,
             "freq"_a = 0.0, "err"_a = 1e-7, "weight"_a = 1.0,
             "blend_tag"_a = "")
        // qn: read returns the full 2*MAXQN flat array; write accepts any
        // Sequence[int] up to length 2*MAXQN and pads the remainder with zeros.
        .def_prop_rw("qn",
            [](const LineRecord &lr) {
                return std::vector<int>(lr.qn.begin(), lr.qn.end());
            },
            [](LineRecord &lr, const std::vector<int> &v) {
                if (v.size() > lr.qn.size())
                    throw nb::value_error("qn list too long (max 2*MAXQN entries)");
                lr.qn.fill(0);
                std::copy(v.begin(), v.end(), lr.qn.begin());
            },
            "Quantum numbers (upper then lower).  Assignment accepts any "
            "Sequence[int] up to 2*MAXQN; shorter lists are padded with zeros.")
        .def_rw("nqn",       &LineRecord::nqn,       "Number of QN pairs used")
        .def_rw("freq",      &LineRecord::freq,      "Observed frequency (MHz)")
        .def_rw("err",       &LineRecord::err,       "Measurement uncertainty (MHz)")
        .def_rw("weight",    &LineRecord::weight,    "Line weight")
        .def_rw("blend_tag", &LineRecord::blend_tag, "Non-empty triggers blend grouping")
        .def("__repr__", [repr_double, repr_string](const LineRecord &lr) {
            std::ostringstream oss;
            // Show only the meaningful slots (2*nqn), or the full array if nqn=0
            int n_show = (lr.nqn > 0) ? std::min(2 * lr.nqn, (int)lr.qn.size())
                                      : (int)lr.qn.size();
            oss << "LineRecord(qn=[";
            for (int i = 0; i < n_show; ++i) {
                if (i) oss << ", ";
                oss << lr.qn[i];
            }
            oss << "], nqn=" << lr.nqn
                << ", freq=" << repr_double(lr.freq)
                << ", err=" << repr_double(lr.err)
                << ", weight=" << repr_double(lr.weight);
            if (!lr.blend_tag.empty())
                oss << ", blend_tag=" << repr_string(lr.blend_tag);
            oss << ")";
            return oss.str();
        });

    nb::class_<DipoleMoment>(m, "DipoleMoment")
        .def("__init__", [](DipoleMoment *self, int64_t id, double value,
                            bool starts_new_component) {
                 new (self) DipoleMoment();
                 self->id = id;
                 self->value = value;
                 self->starts_new_component = starts_new_component;
             },
             "id"_a = 0, "value"_a = 0.0, "starts_new_component"_a = false)
        .def_rw("id",                   &DipoleMoment::id,                   "Dipole identifier (decimal)")
        .def_rw("value",                &DipoleMoment::value,                "Dipole moment (Debye)")
        .def_rw("starts_new_component", &DipoleMoment::starts_new_component, "True → NEGBCD flag (new component group)")
        .def("__repr__", [repr_double](const DipoleMoment &d) {
            std::ostringstream oss;
            oss << "DipoleMoment(id=" << d.id
                << ", value=" << repr_double(d.value)
                << ", starts_new_component=" << (d.starts_new_component ? "True" : "False") << ")";
            return oss.str();
        });

    nb::class_<VibState>(m, "VibState")
        .def("__init__", [](VibState *self, int index, int knmin, int knmax,
                            int stat_weight_axis, int iwtpl, int iwtmn,
                            double vsym, int esym_weight,
                            std::vector<int> spin_degeneracies,
                            bool symmetric_rotor_quanta) {
                 new (self) VibState();
                 self->index = index;
                 self->knmin = knmin;
                 self->knmax = knmax;
                 self->stat_weight_axis = stat_weight_axis;
                 self->iwtpl = iwtpl;
                 self->iwtmn = iwtmn;
                 self->vsym = vsym;
                 self->esym_weight = esym_weight;
                 self->spin_degeneracies = std::move(spin_degeneracies);
                 self->symmetric_rotor_quanta = symmetric_rotor_quanta;
             },
             "index"_a = 0, "knmin"_a = 0, "knmax"_a = 359,
             "stat_weight_axis"_a = 1, "iwtpl"_a = 1, "iwtmn"_a = 1,
             "vsym"_a = 0.0, "esym_weight"_a = 99,
             "spin_degeneracies"_a = std::vector<int>{},
             "symmetric_rotor_quanta"_a = false)
        .def_rw("index",                  &VibState::index,                  "0-based vibrational state index")
        .def_rw("knmin",                  &VibState::knmin,                  "Minimum K quantum number (set both to 0 for linear molecules)")
        .def_rw("knmax",                  &VibState::knmax,                  "Maximum K quantum number (set both to 0 for linear molecules)")
        .def_rw("stat_weight_axis",       &VibState::stat_weight_axis,       "Axis for statistical weight: 1=a, 2=b, 3=c, 4=A(C2)... Negative → I_tot spin basis")
        .def_rw("iwtpl",                  &VibState::iwtpl,                  "Statistical weight for even (plus) parity states")
        .def_rw("iwtmn",                  &VibState::iwtmn,                  "Statistical weight for odd (minus) parity states")
        .def_rw("vsym",                   &VibState::vsym,                   "Vibrational symmetry (only meaningful on last VibState; builder injects -1.0 sentinel for earlier states)")
        .def_rw("esym_weight",            &VibState::esym_weight,            "E-symmetry weight (EWT). Default 99 = ignored. Negative triggers EWTFAC=1000 mode.")
        .def_rw("spin_degeneracies",      &VibState::spin_degeneracies,      "Spin degeneracy (2I+1 or 2S+1) for each spin species. Empty = no spin species.")
        .def_rw("symmetric_rotor_quanta", &VibState::symmetric_rotor_quanta, "True for linear molecules and open-shell systems (selects symmetric-top K quantum numbers)")
        .def("__repr__", [repr_double](const VibState &v) {
            std::ostringstream oss;
            oss << "VibState(index=" << v.index
                << ", knmin=" << v.knmin << ", knmax=" << v.knmax
                << ", stat_weight_axis=" << v.stat_weight_axis
                << ", iwtpl=" << v.iwtpl << ", iwtmn=" << v.iwtmn
                << ", vsym=" << repr_double(v.vsym)
                << ", esym_weight=" << v.esym_weight
                << ", spin_degeneracies=[";
            for (size_t i = 0; i < v.spin_degeneracies.size(); ++i) {
                if (i) oss << ", ";
                oss << v.spin_degeneracies[i];
            }
            oss << "], symmetric_rotor_quanta=" << (v.symmetric_rotor_quanta ? "True" : "False") << ")";
            return oss.str();
        });

    nb::class_<SpinvOptions>(m, "SpinvOptions")
        .def("__init__", [](SpinvOptions *self, int inclusion_flags, int diag_order,
                            int phase_flags, bool oblate, std::string nam_file,
                            std::vector<VibState> vibs) {
                 new (self) SpinvOptions();
                 self->inclusion_flags = inclusion_flags;
                 self->diag_order = diag_order;
                 self->phase_flags = phase_flags;
                 self->oblate = oblate;
                 self->nam_file = std::move(nam_file);
                 self->vibs = std::move(vibs);
             },
             "inclusion_flags"_a = 0, "diag_order"_a = 0, "phase_flags"_a = 0,
             "oblate"_a = false, "nam_file"_a = "",
             "vibs"_a = std::vector<VibState>{})
        .def_rw("inclusion_flags", &SpinvOptions::inclusion_flags, "Binary flags for inter-state couplings: bit 0 = no ΔN≠0, bit 1 = no ΔJ, etc. 0 = include all.")
        .def_rw("diag_order",      &SpinvOptions::diag_order,      "Eigenvalue ordering: 0=energy, 1=full projection, 2=energy/Wang, 3=τ=Ka-Kc, 4=<Kz²>, 5=diag")
        .def_rw("phase_flags",     &SpinvOptions::phase_flags,     "Packed phase flags: PHASE + 10*NEWLZ + 20*NOFC + 40*G12")
        .def_rw("oblate",          &SpinvOptions::oblate,          "Oblate rotor: z=c, y=b, x=a (default False = prolate: z=a)")
        .def_rw("nam_file",        &SpinvOptions::nam_file,        "Parameter-label file path (empty = engine default 'sping.nam')")
        .def_rw("vibs",            &SpinvOptions::vibs,            "Vibrational state list (at least one entry required)")
        .def("__repr__", [repr_string](const SpinvOptions &s) {
            std::ostringstream oss;
            oss << "SpinvOptions(inclusion_flags=" << s.inclusion_flags
                << ", diag_order=" << s.diag_order
                << ", phase_flags=" << s.phase_flags
                << ", oblate=" << (s.oblate ? "True" : "False")
                << ", nam_file=" << repr_string(s.nam_file)
                << ", vibs=[" << s.vibs.size() << " VibState"
                << (s.vibs.size() == 1 ? "" : "s") << "])";
            return oss.str();
        });

    nb::class_<DpiOptions>(m, "DpiOptions")
        .def("__init__", [](DpiOptions *self, int isdgn, int nvib) {
                 new (self) DpiOptions();
                 self->isdgn = isdgn;
                 self->nvib = nvib;
             },
             "isdgn"_a = 1, "nvib"_a = 1)
        .def_rw("isdgn", &DpiOptions::isdgn, "Spin degeneracy")
        .def_rw("nvib",  &DpiOptions::nvib,  "Number of vibrational states")
        .def("__repr__", [](const DpiOptions &d) {
            std::ostringstream oss;
            oss << "DpiOptions(isdgn=" << d.isdgn << ", nvib=" << d.nvib << ")";
            return oss.str();
        });

    // EngineOptions: kind is inferred from which sub-options are provided.
    //   EngineOptions()                        → kind=Spinv
    //   EngineOptions(spinv=...)               → kind=Spinv
    //   EngineOptions(dpi=...)                 → kind=Dpi
    //   EngineOptions(kind=..., spinv|dpi=...) → explicit override
    nb::class_<EngineOptions>(m, "EngineOptions")
        .def("__init__", [](EngineOptions *self,
                            nb::object kind_obj,
                            nb::object spinv_obj,
                            nb::object dpi_obj) {
                 new (self) EngineOptions();
                 bool has_spinv = !spinv_obj.is_none();
                 bool has_dpi   = !dpi_obj.is_none();
                 if (!kind_obj.is_none()) {
                     self->kind = nb::cast<EngineOptions::Kind>(kind_obj);
                 } else if (has_dpi && !has_spinv) {
                     self->kind = EngineOptions::Kind::Dpi;
                 } else {
                     self->kind = EngineOptions::Kind::Spinv;
                 }
                 if (has_spinv) self->spinv = nb::cast<SpinvOptions>(spinv_obj);
                 if (has_dpi)   self->dpi   = nb::cast<DpiOptions>(dpi_obj);
             },
             "kind"_a = nb::none(), "spinv"_a = nb::none(), "dpi"_a = nb::none())
        .def_rw("kind",  &EngineOptions::kind,  "Which engine to use (EngineKind.Spinv or EngineKind.Dpi)")
        .def_rw("spinv", &EngineOptions::spinv, "SPINV engine options (valid when kind=Spinv)")
        .def_rw("dpi",   &EngineOptions::dpi,   "DPI engine options (valid when kind=Dpi)")
        .def("__repr__", [](const EngineOptions &e) {
            std::ostringstream oss;
            oss << "EngineOptions(kind=EngineKind."
                << (e.kind == EngineOptions::Kind::Spinv ? "Spinv" : "Dpi") << ")";
            return oss.str();
        });

    nb::class_<CatControl>(m, "CatControl")
        .def("__init__", [](CatControl *self, int iflg, long itag, double qrot,
                            int inblk, int lblk, double thrsh, double thrsh1,
                            double fqmax, double tmq, int maxv) {
                 new (self) CatControl();
                 self->iflg = iflg;
                 self->itag = itag;
                 self->qrot = qrot;
                 self->inblk = inblk;
                 self->lblk = lblk;
                 self->thrsh = thrsh;
                 self->thrsh1 = thrsh1;
                 self->fqmax = fqmax;
                 self->tmq = tmq;
                 self->maxv = maxv;
             },
             "iflg"_a = 0, "itag"_a = 999L, "qrot"_a = 1000.0,
             "inblk"_a = 0, "lblk"_a = 0, "thrsh"_a = -100.0, "thrsh1"_a = -100.0,
             "fqmax"_a = 9999.99, "tmq"_a = 300.0, "maxv"_a = -1)
        .def_rw("iflg",   &CatControl::iflg,   "Packed flags: IRFLG*1000 + OUTFLG*100 + STRFLG*10 + EGYFLG")
        // Convenience accessors for the iflg packed fields:
        .def_prop_rw("wavenumbers",
            [](const CatControl &c) { return (c.iflg / 1000) % 10 != 0; },
            [](CatControl &c, bool v) {
                c.iflg = (c.iflg % 1000) + (v ? 1000 : 0);
            },
            "True if parameters are in wavenumbers (sets IRFLG)")
        .def_prop_rw("output_strengths",
            [](const CatControl &c) { return (c.iflg / 10) % 10; },
            [](CatControl &c, int v) {
                if (v < 0 || v > 9)
                    throw nb::value_error("output_strengths must be 0-9");
                c.iflg = (c.iflg / 100) * 100 + v * 10 + (c.iflg % 10);
            },
            "STRFLG: 0=off, 1=enable .str, 2=label separate dipole contributions")
        .def_prop_rw("output_energies",
            [](const CatControl &c) { return c.iflg % 10; },
            [](CatControl &c, int v) {
                if (v < 0 || v > 9)
                    throw nb::value_error("output_energies must be 0-9");
                c.iflg = (c.iflg / 10) * 10 + v;
            },
            "EGYFLG: 0=off, 1=energies, 2=+derivatives, 3=+eigenvectors, 5=dump Hamiltonian")
        .def_rw("itag",   &CatControl::itag,   "Species tag")
        .def_rw("qrot",   &CatControl::qrot,   "Partition function scaling")
        .def_rw("inblk",  &CatControl::inblk,  "Beginning F quantum (rounded up)")
        .def_rw("lblk",   &CatControl::lblk,   "Ending F quantum")
        .def_rw("thrsh",  &CatControl::thrsh,  "Intensity threshold (log10 scale)")
        .def_rw("thrsh1", &CatControl::thrsh1, "Secondary intensity threshold")
        .def_rw("fqmax",  &CatControl::fqmax,  "Maximum frequency (GHz)")
        .def_rw("tmq",    &CatControl::tmq,    "Temperature for intensity (K)")
        .def_rw("maxv",   &CatControl::maxv,   "Maximum v quantum number (-1 = unlimited)")
        .def("__repr__", [repr_double](const CatControl &c) {
            std::ostringstream oss;
            oss << "CatControl(iflg=" << c.iflg
                << ", itag=" << c.itag
                << ", qrot=" << repr_double(c.qrot)
                << ", fqmax=" << repr_double(c.fqmax)
                << ", tmq=" << repr_double(c.tmq)
                << ", maxv=" << c.maxv << ")";
            return oss.str();
        });

    nb::class_<FitInput>(m, "FitInput")
        .def("__init__", [](FitInput *self,
                            std::string title, int n_iterations,
                            double marquardt_param, double max_obs_calc_err,
                            double param_err_scale, double freq_scale,
                            size_t max_lines, bool extended_qn, int nxpar,
                            EngineOptions engine_options,
                            std::vector<Parameter> parameters,
                            std::vector<double> variance,
                            std::vector<LineRecord> lines) {
                 new (self) FitInput();
                 self->title = std::move(title);
                 self->n_iterations = n_iterations;
                 self->marquardt_param = marquardt_param;
                 self->max_obs_calc_err = max_obs_calc_err;
                 self->param_err_scale = param_err_scale;
                 self->freq_scale = freq_scale;
                 self->max_lines = max_lines;
                 self->extended_qn = extended_qn;
                 self->nxpar = nxpar;
                 self->engine_options = std::move(engine_options);
                 self->parameters = std::move(parameters);
                 self->variance = std::move(variance);
                 self->lines = std::move(lines);
             },
             "title"_a = "", "n_iterations"_a = 1,
             "marquardt_param"_a = 0.0, "max_obs_calc_err"_a = 1e6,
             "param_err_scale"_a = 1.0, "freq_scale"_a = 1.0,
             "max_lines"_a = (size_t)32767, "extended_qn"_a = false, "nxpar"_a = 0,
             "engine_options"_a = EngineOptions{},
             "parameters"_a = std::vector<Parameter>{},
             "variance"_a = std::vector<double>{},
             "lines"_a = std::vector<LineRecord>{})
        .def_rw("title",            &FitInput::title)
        .def_rw("n_iterations",     &FitInput::n_iterations,     "Number of fit iterations")
        .def_rw("marquardt_param",  &FitInput::marquardt_param,  "Marquardt damping parameter")
        .def_rw("max_obs_calc_err", &FitInput::max_obs_calc_err, "Maximum obs-calc/err to include")
        .def_rw("param_err_scale",  &FitInput::param_err_scale,  "Parameter error scale factor (parfac)")
        .def_rw("freq_scale",       &FitInput::freq_scale,       "Frequency scale factor (fqfacq)")
        .def_rw("max_lines",        &FitInput::max_lines,        "Maximum number of lines")
        .def_rw("extended_qn",      &FitInput::extended_qn,      "Use MAXQN QNs instead of MAXCAT")
        .def_rw("nxpar",            &FitInput::nxpar,            "Number of parameters excluded if index < 0")
        .def_rw("engine_options",   &FitInput::engine_options)
        .def_rw("parameters",       &FitInput::parameters)
        .def_rw("variance",         &FitInput::variance,         "Packed upper-triangular variance matrix")
        .def_rw("lines",            &FitInput::lines)
        .def("__repr__", [repr_string](const FitInput &fi) {
            std::ostringstream oss;
            oss << "FitInput(title=" << repr_string(fi.title)
                << ", parameters=" << fi.parameters.size()
                << ", lines=" << fi.lines.size()
                << ", kind=" << (fi.engine_options.kind == EngineOptions::Kind::Spinv ? "Spinv" : "Dpi")
                << ")";
            return oss.str();
        });

    nb::class_<CatInput>(m, "CatInput")
        .def("__init__", [](CatInput *self,
                            std::string title, CatControl control,
                            std::vector<DipoleMoment> dipoles,
                            EngineOptions engine_options,
                            std::vector<Parameter> parameters,
                            std::vector<double> variance) {
                 new (self) CatInput();
                 self->title = std::move(title);
                 self->control = std::move(control);
                 self->dipoles = std::move(dipoles);
                 self->engine_options = std::move(engine_options);
                 self->parameters = std::move(parameters);
                 self->variance = std::move(variance);
             },
             "title"_a = "", "control"_a = CatControl{},
             "dipoles"_a = std::vector<DipoleMoment>{},
             "engine_options"_a = EngineOptions{},
             "parameters"_a = std::vector<Parameter>{},
             "variance"_a = std::vector<double>{})
        .def_rw("title",          &CatInput::title)
        .def_rw("control",        &CatInput::control)
        .def_rw("dipoles",        &CatInput::dipoles)
        .def_rw("engine_options", &CatInput::engine_options)
        .def_rw("parameters",     &CatInput::parameters)
        .def_rw("variance",       &CatInput::variance, "Packed upper-triangular variance matrix")
        .def("__repr__", [repr_string](const CatInput &ci) {
            std::ostringstream oss;
            oss << "CatInput(title=" << repr_string(ci.title)
                << ", parameters=" << ci.parameters.size()
                << ", dipoles=" << ci.dipoles.size()
                << ", kind=" << (ci.engine_options.kind == EngineOptions::Kind::Spinv ? "Spinv" : "Dpi")
                << ")";
            return oss.str();
        });

    // ---- Typed-struct parsers ----
    m.def("parse_fit_files", &parse_fit_files,
          "par_file"_a, "lin_file"_a, "kind"_a = EngineOptions::Kind::Spinv,
          "Parse .par and .lin files into a FitInput struct (no engine configured).");

    m.def("parse_cat_files", &parse_cat_files,
          "var_file"_a, "int_file"_a, "kind"_a = EngineOptions::Kind::Spinv,
          "Parse .var and .int files into a CatInput struct (no engine configured).");
}
