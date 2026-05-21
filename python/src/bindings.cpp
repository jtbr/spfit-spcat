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
        .def_ro("par",     &CalFitOutput::par,     "Fitted parameter values")
        .def_ro("erpar",   &CalFitOutput::erpar,   "Estimated parameter errors")
        .def_ro("xsqbest", &CalFitOutput::xsqbest, "Best RMS (obs-calc)/err after fitting")
        .def_ro("itr",     &CalFitOutput::itr,     "Number of iterations performed");

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
        }, "Partition function Q at each temperature");

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

    nb::enum_<EngineOptions::Kind>(m, "EngineKind")
        .value("Spinv", EngineOptions::Kind::Spinv)
        .value("Dpi",   EngineOptions::Kind::Dpi);

    nb::class_<Parameter>(m, "Parameter")
        .def(nb::init<>())
        .def_rw("id",             &Parameter::id,             "Decimal parameter ID")
        .def_rw("value",          &Parameter::value,          "Absolute value (MHz or cm-1)")
        .def_rw("a_priori_error", &Parameter::a_priori_error, "1-sigma uncertainty; <=1e-37 fixes parameter")
        .def_rw("fixed",          &Parameter::fixed,          "Excluded from fit if true")
        .def_rw("label",          &Parameter::label,          "Optional label (truncated to 10 chars)");

    nb::class_<LineRecord>(m, "LineRecord")
        .def(nb::init<>())
        .def_rw("qn",        &LineRecord::qn,        "Upper then lower QNs (2*MAXQN entries)")
        .def_rw("nqn",       &LineRecord::nqn,       "Number of QN pairs used")
        .def_rw("freq",      &LineRecord::freq,      "Observed frequency (MHz)")
        .def_rw("err",       &LineRecord::err,       "Measurement uncertainty (MHz)")
        .def_rw("weight",    &LineRecord::weight,    "Line weight")
        .def_rw("blend_tag", &LineRecord::blend_tag, "Non-empty triggers blend grouping");

    nb::class_<DipoleMoment>(m, "DipoleMoment")
        .def(nb::init<>())
        .def_rw("id",                   &DipoleMoment::id,                   "Dipole identifier (decimal)")
        .def_rw("value",                &DipoleMoment::value,                "Dipole moment (Debye)")
        .def_rw("starts_new_component", &DipoleMoment::starts_new_component, "True → NEGBCD flag (new component group)");

    nb::class_<VibState>(m, "VibState")
        .def(nb::init<>())
        .def_rw("index",         &VibState::index,         "0-based vibrational state index")
        .def_rw("knmin",         &VibState::knmin,         "Minimum K quantum number")
        .def_rw("knmax",         &VibState::knmax,         "Maximum K (0 → linear molecule)")
        .def_rw("iax",           &VibState::iax,           "Axis / symmetry selector (0-11)")
        .def_rw("iwtpl",         &VibState::iwtpl,         "Statistical weight (plus parity)")
        .def_rw("iwtmn",         &VibState::iwtmn,         "Statistical weight (minus parity)")
        .def_rw("vsym",          &VibState::vsym,          "Vibrational symmetry parameter")
        .def_rw("ewt0",          &VibState::ewt0,          "Energy/weight flag")
        .def_rw("nuclear_spins", &VibState::nuclear_spins, "2*I values for each nucleus")
        .def_rw("negbcd_spin",   &VibState::negbcd_spin,   "NEGBCD flag on spin field (required for linear molecules)");

    nb::class_<SpinvOptions>(m, "SpinvOptions")
        .def(nb::init<>())
        .def_rw("ixz",         &SpinvOptions::ixz,         "x-z plane selection")
        .def_rw("idiag",       &SpinvOptions::idiag,       "Diagonalisation option")
        .def_rw("phase_flags", &SpinvOptions::phase_flags, "Packed phase flags (stdphase/newlz/nofc/g12)")
        .def_rw("oblate",      &SpinvOptions::oblate,      "Global oblate flag")
        .def_rw("nam_file",    &SpinvOptions::nam_file,    "Parameter-name file path (empty = engine default)")
        .def_rw("vibs",        &SpinvOptions::vibs,        "Vibrational state list (at least one required)");

    nb::class_<DpiOptions>(m, "DpiOptions")
        .def(nb::init<>())
        .def_rw("isdgn", &DpiOptions::isdgn, "Spin degeneracy")
        .def_rw("nvib",  &DpiOptions::nvib,  "Number of vibrational states");

    nb::class_<EngineOptions>(m, "EngineOptions")
        .def(nb::init<>())
        .def_rw("kind",  &EngineOptions::kind,  "Which engine to use (EngineKind.Spinv or EngineKind.Dpi)")
        .def_rw("spinv", &EngineOptions::spinv, "SPINV engine options (valid when kind=Spinv)")
        .def_rw("dpi",   &EngineOptions::dpi,   "DPI engine options (valid when kind=Dpi)");

    nb::class_<CatControl>(m, "CatControl")
        .def(nb::init<>())
        .def_rw("iflg",   &CatControl::iflg)
        .def_rw("itag",   &CatControl::itag,   "Species tag")
        .def_rw("qrot",   &CatControl::qrot,   "Partition function scaling")
        .def_rw("inblk",  &CatControl::inblk)
        .def_rw("lblk",   &CatControl::lblk)
        .def_rw("thrsh",  &CatControl::thrsh,  "Intensity threshold (log10 scale)")
        .def_rw("thrsh1", &CatControl::thrsh1)
        .def_rw("fqmax",  &CatControl::fqmax,  "Maximum frequency (GHz)")
        .def_rw("tmq",    &CatControl::tmq,    "Temperature for intensity (K)")
        .def_rw("maxv",   &CatControl::maxv,   "Maximum v quantum number (-1 = unlimited)");

    nb::class_<FitInput>(m, "FitInput")
        .def(nb::init<>())
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
        .def_rw("lines",            &FitInput::lines);

    nb::class_<CatInput>(m, "CatInput")
        .def(nb::init<>())
        .def_rw("title",          &CatInput::title)
        .def_rw("control",        &CatInput::control)
        .def_rw("dipoles",        &CatInput::dipoles)
        .def_rw("engine_options", &CatInput::engine_options)
        .def_rw("parameters",     &CatInput::parameters)
        .def_rw("variance",       &CatInput::variance, "Packed upper-triangular variance matrix");

    // ---- Typed-struct parsers ----
    m.def("parse_fit_files", &parse_fit_files,
          "par_file"_a, "lin_file"_a, "kind"_a = EngineOptions::Kind::Spinv,
          "Parse .par and .lin files into a FitInput struct (no engine configured).");

    m.def("parse_cat_files", &parse_cat_files,
          "var_file"_a, "int_file"_a, "kind"_a = EngineOptions::Kind::Spinv,
          "Parse .var and .int files into a CatInput struct (no engine configured).");
}
