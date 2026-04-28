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
        "Read SPFIT input files and run the fitter.  Each instance may call run() once.")
        .def(nb::init<const std::string &, const std::string &, const std::string &>(),
             "par_file"_a, "lin_file"_a, "engine"_a = "spinv",
             "Open par_file and lin_file; configure the engine via setopt.")
        .def("run", &FitSession::run,
             "Run the fitting and return a CalFitOutput.  May only be called once.");

    nb::class_<CatSession>(m, "CatSession",
        "Read SPCAT input files and generate the catalog.  Each instance may call run() once.")
        .def(nb::init<const std::string &, const std::string &, const std::string &>(),
             "int_file"_a, "var_file"_a, "engine"_a = "spinv",
             "Open int_file and var_file; configure the engine via setopt.")
        .def("run", &CatSession::run,
             "Generate catalog and return a CalCatOutput.  May only be called once.");

    // ---- High-level convenience functions ----
    m.def("fit_from_files", &fit_from_files,
          "par_file"_a, "lin_file"_a, "engine"_a = "spinv",
          "Read SPFIT input files, run the fitter, and return CalFitOutput.");

    m.def("cat_from_files", &cat_from_files,
          "int_file"_a, "var_file"_a, "engine"_a = "spinv",
          "Read SPCAT input files, generate the catalog, and return CalCatOutput.");
}
