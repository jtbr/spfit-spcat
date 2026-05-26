/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */

/**************************************************************************/
/*                                                                        */
/*   SPFIT — Spectroscopic Parameter Fitting Program                      */
/*   THIS IS A GENERALIZED LINE FITTING PROGRAM                           */
/*   IT FITS LINES TO PARAMETERS IN A MODEL HAMILTONIAN                   */
/*   BLENDED LINES ARE TREATED SPECIALLY:                                 */
/*     IF THE EXPTL.FREQ. IS THE SAME TO 1 HZ THEN ONLY THE INVERSE ERROR */
/*             AVERAGED FREQUENCIES ARE USED IN THE FIT FOR THE BLEND     */
/*                                                                        */
/**************************************************************************/
/* "trust region" Marquardt fitting is described in John. E. Dennis and   */
/* Robert B. Schnabel, Numerical Methods for Unsconstrained Optimization  */
/* and Non-linear Equations, Prentice-Hall, 1983.                         */

#include <stdio.h>
#include <string.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#include "splib/lsqfit.h"
#include "spfit/subfit.h"
#include "engine/SpinvEngine.hpp"
#include "engine/DpiEngine.hpp"
#include "CalFit.hpp"
#include "CalFitIO.hpp"
#include "api/toml_io.hpp"
#include "api/builders.hpp"
#include "common/CalError.hpp"
#include "common/file_helpers.hpp"
#include "common/SigintFlag.hpp"
#include "api/legacy_parser.hpp"

/**
 * @brief Main function for SPFIT - Spectroscopic Parameter Fitting program
 *
 * This program fits spectroscopic parameters to experimental line frequencies.
 * It uses the CalFit class to perform the fitting process, delegating the core
 * logic to a modular implementation.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @return int Exit status (0 for success)
 */
int main(int argc, char *argv[])
{
#define NFILE 5
  static const char *ext[NFILE] = {"par", "lin", "fit", "bak", "var"};
  enum efile { epar, elin, efit, ebak, evar };
  char *fname[NFILE + 1];
  std::string engineType = "spinv";
  bool toml_out = false;
  // Strip leading option flags in any order: --spinv | --dpi | --toml-out
  for (bool stripped = true; stripped && argc > 1; ) {
    stripped = false;
    if      (strcasecmp(argv[1], "--dpi")      == 0) { engineType = "dpi"; stripped = true; }
    else if (strcasecmp(argv[1], "--spinv")    == 0) {                      stripped = true; }
    else if (strcasecmp(argv[1], "--toml-out") == 0) { toml_out = true;     stripped = true; }
    if (stripped) { char *n = *argv; argv++; *argv = n; argc--; }
  }

  Logger &logger = Logger::defaultLogger();
  logger.info("Using %s engine.", engineType.c_str());

  file_helpers::parse_file_args(argc, argv, NFILE, fname, ext);

  // ── TOML mode: mol.toml present ──────────────────────────────────────────
  //
  // Derive base from fname[epar] ("base.par" → "base").
  // If "base.toml" exists, use the TOML path instead of the legacy .par/.lin path.

  char *base = file_helpers::base_name(fname[epar]);
  std::string toml_in  = base ? (std::string(base) + ".toml")     : "";
  std::string var_toml = base ? (std::string(base) + ".var.toml") : "";
  free(base);

  bool toml_mode = !toml_in.empty() && file_helpers::file_exists(toml_in.c_str());

  if (toml_mode) {
    // ── TOML path ─────────────────────────────────────────────────────────
    FILE *lufit = fopen(fname[efit], "w");
    if (!lufit) {
      fprintf(stderr, "spfit: cannot open %s\n", fname[efit]);
      return EXIT_FAILURE;
    }

    try {
      std::unique_ptr<CalculationEngine> eng =
          (engineType == "dpi") ? std::unique_ptr<CalculationEngine>(std::make_unique<DpiEngine>())
                                : std::unique_ptr<CalculationEngine>(std::make_unique<SpinvEngine>());

      FitInput fi = load_fit_input_toml(toml_in);
      // Override engine kind from CLI flag if given
      if (engineType == "dpi")   fi.engine_options.kind = EngineOptions::Kind::Dpi;
      if (engineType == "spinv") fi.engine_options.kind = EngineOptions::Kind::Spinv;

      CalFitInput ci = build_fit_input(fi, *eng);
      CalFit calFit(eng, lufit, logger);
      CalFitOutput output;
      SigintFlag sigint_guard;
      calFit.run(ci, output);
      fclose(lufit);

      save_fit_output_toml(output, fi, var_toml);
      logger.info("FIT COMPLETE (TOML mode) → %s", var_toml.c_str());
    } catch (const CalError &e) {
      fprintf(stderr, "spfit error: %s\n", e.what());
      fclose(lufit);
      return EXIT_FAILURE;
    }
    return 0;
  }

  // ── Legacy path (.par / .lin) ─────────────────────────────────────────────

  std::unique_ptr<CalculationEngine> calc_engine =
      (engineType == "dpi") ? std::unique_ptr<CalculationEngine>(std::make_unique<DpiEngine>())
                            : std::unique_ptr<CalculationEngine>(std::make_unique<SpinvEngine>());

  // Backup .par file to .bak file
  // Original: filbak(fname[epar], fname[ebak]);
  // fname[ebak] is "base.bak"
  // This means "base.par" is copied to "base.bak"
  // The program then reads from "base.par" (which is now fname[ebak] in original logic if lubak was fname[ebak])
  // And writes the new .par to fname[epar].
  // Let's stick to: read from original fname[epar], backup fname[epar] to fname[ebak].
  // CalFitIO::writeOutput will write to a new fname[epar].

  if (filbak(fname[epar], fname[ebak]))
  { // filbak returns 0 on success
    logger.warn("Failed to create backup file %s from %s.", fname[ebak], fname[epar]);
    // Or exit if backup is critical
  } else {
    logger.info("Successfully created backup file %s from %s.", fname[ebak], fname[epar]);
  }

  // Open the main .fit output stream
  FILE *lufit_stream = fopen(fname[efit], "w");
  if (!lufit_stream) {
    logger.error("Unable to open fit output file '%s'.", fname[efit]);
    return EXIT_FAILURE;
  }

  // Prepare input and output structures
  CalFitInput input;
  CalFitOutput output;

  // Read input data using CalFitIO
  // We read from the original .par file (fname[epar])
  if (!CalFitIO::readInput(fname[ebak], fname[elin], input, calc_engine, lufit_stream)) {
    logger.error("Failed to read input files.");
    fclose(lufit_stream);
    return EXIT_FAILURE;
  }

  // Create CalFit instance with the selected engine and lufit stream
  // (calFit takes ownership of calc_engine from here)
  CalFit calFit(calc_engine, lufit_stream, logger);

  // Run the fitting process (guard active for duration of run)
  SigintFlag sigint_guard;
  try {
    calFit.run(input, output);
  } catch (const CalError &e) {
    logger.error("Fitting process failed: %s", e.what());
    fclose(lufit_stream);
    return EXIT_FAILURE;
  }

  fclose(lufit_stream);

  // Write output data using CalFitIO
  // CalFitIO::writeOutput writes to:
  // - fname[epar] (new parameter file)
  // - fname[evar] (new variance file)
  // fname[ebak] (the backup of the original .par) is kept.
  if (!CalFitIO::writeOutput(fname[epar], fname[ebak], fname[evar], output, input))
  {
    logger.error("Failed to write output files.");
    return EXIT_FAILURE;
  }

  if (toml_out && !var_toml.empty()) {
    try {
      FitInput fi = parse_fit_files(fname[ebak], fname[elin]);
      if (engineType == "dpi") fi.engine_options.kind = EngineOptions::Kind::Dpi;
      save_fit_output_toml(output, fi, var_toml);
      logger.info("--toml-out: wrote %s", var_toml.c_str());
    } catch (const CalError &e) {
      logger.warn("--toml-out: failed to write %s: %s", var_toml.c_str(), e.what());
    }
  }

  logger.info("FIT COMPLETE");
  return 0;
}