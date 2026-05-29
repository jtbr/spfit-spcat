/***********************************************************************/
/*   Copyright (C) 1989, California Institute of Technology            */
/*   All rights reserved.  U. S. Government Sponsorship under          */
/*   NASA Contract NAS7-918 is acknowledged.                           */

/*   Herbert M. Pickett, 20 Mar 1989                                   */
/*   Revised version in c, 22 March 1999                               */
/*   Revised version in C++, 2026                                      */

/* SPCAT — Spectroscopic Catalog Generator                             */

/* THIS IS A GENERALIZED INTENSITY AND FREQUENCY CALCULATOR            */
/* THE HAMILTONIAN IS ASSUMED TO BE ORGANIZED INTO BLOCKS              */
/*     SUCH THAT ADJACENT SETS OF BLOCKS ARE CONNECTED BY TRANSITIONS  */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#  define strcasecmp _stricmp
#endif
#include "engine/SpinvEngine.hpp"
#include "engine/DpiEngine.hpp"
#include "CalCat.hpp"
#include "CalCatIO.hpp"
#include "splib/calpgm_types.h"
#include "splib/ulib.h"
#include "spcat/OutputSink.hpp"
#include "api/toml_io.hpp"
#include "api/builders.hpp"
#include "common/CalError.hpp"
#include "common/file_helpers.hpp"
#include "common/SigintFlag.hpp"

int main(int argc, char *argv[])
{
#define NFILE 6
  static const char *ext[NFILE] =
      {"int", "var", "out", "cat", "str", "egy"};
  enum efile { eint, evar, eout, ecat, estr, eegy };
  char *fname[NFILE + 1];

  /* Choose the calculation engine (SPINV by default, or DPI) */
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

  file_helpers::parse_file_args(argc, argv, NFILE, fname, ext);

  // ── TOML mode: mol.fitted.toml + mol.dipoles.toml present ───────────────
  //
  // Derive base from fname[evar] ("base.var" → "base").
  // If "base.fitted.toml" exists, use the TOML path.

  char *base = file_helpers::base_name(fname[evar]);
  std::string base_str  = base ? std::string(base) : "";
  std::string var_toml  = base_str.empty() ? "" : (base_str + ".fitted.toml");
  std::string int_toml  = base_str.empty() ? "" : (base_str + ".dipoles.toml");
  std::string cat_toml  = base_str.empty() ? "" : (base_str + ".catalog.toml");
  free(base);

  // Also accept dipoles.toml derived from fname[eint]
  if (!int_toml.empty() && !file_helpers::file_exists(int_toml.c_str())) {
    char *ibase = file_helpers::base_name(fname[eint]);
    if (ibase) { int_toml = std::string(ibase) + ".dipoles.toml"; free(ibase); }
  }

  bool toml_mode = !var_toml.empty() && file_helpers::file_exists(var_toml.c_str())
                && !int_toml.empty() && file_helpers::file_exists(int_toml.c_str());

  fprintf(stderr, "Using %s engine, %s mode.\n", engineType.c_str(), toml_mode ? "TOML" : "legacy");
  if (toml_mode) {
    if (file_helpers::file_exists(fname[evar]))
      fprintf(stderr, "WARNING: TOML mode active but legacy file %s also present (ignored).\n", fname[evar]);
    if (file_helpers::file_exists(fname[eint]))
      fprintf(stderr, "WARNING: TOML mode active but legacy file %s also present (ignored).\n", fname[eint]);
  }

  if (toml_mode) {
    // ── TOML path ─────────────────────────────────────────────────────────
    try {
      CatInput ci = load_cat_input_toml(var_toml, int_toml);
      if (engineType == "dpi")   ci.engine_options.kind = EngineOptions::Kind::Dpi;
      if (engineType == "spinv") ci.engine_options.kind = EngineOptions::Kind::Spinv;

      std::unique_ptr<CalculationEngine> eng =
          (ci.engine_options.kind == EngineOptions::Kind::Dpi)
          ? std::unique_ptr<CalculationEngine>(std::make_unique<DpiEngine>())
          : std::unique_ptr<CalculationEngine>(std::make_unique<SpinvEngine>());

      CalCatInput cci = build_cat_input(ci, *eng);

      FILE *luout_f = file_helpers::open_output(fname[eout], "w");
      CalCatIO::write_cat_preamble(luout_f, ci, cci);

      // Route .egy and .str to separate files if flags request it (matching legacy path)
      int iflg_tmp = cci.iflg;
      if (iflg_tmp >= 1000) iflg_tmp %= 1000;
      iflg_tmp %= 100;
      FILE *luegy_f = luout_f;
      FILE *lustr_f = luout_f;
      std::string egy_path = base_str.empty() ? "" : (base_str + ".egy");
      std::string str_path = base_str.empty() ? "" : (base_str + ".str");
      if (iflg_tmp >= 10 && !str_path.empty())
        lustr_f = file_helpers::open_output(str_path.c_str(), "w");
      if ((iflg_tmp % 10) != 0 && !egy_path.empty())
        luegy_f = file_helpers::open_output(egy_path.c_str(), "w");

      FileSink luout_sink(luout_f);
      MemorySink lucat_sink;
      FileSink   luegy_sink(luegy_f);
      FileSink   lustr_sink(lustr_f);

      CalCatOutput output;
      {
        SigintFlag sigint_guard;
        CalCat calCat(eng, &luout_sink, &lucat_sink, &luegy_sink, &lustr_sink);
        calCat.run(cci, output);
      }

      output.cat_lines = lucat_sink.drain_lines();
      output.sort_cat_lines();
      if (luegy_f != luout_f) fclose(luegy_f);
      if (lustr_f != luout_f) fclose(lustr_f);
      fclose(luout_f);

      save_cat_output_toml(output, cat_toml);
      fprintf(stderr, "spcat: wrote %s\n", cat_toml.c_str());
    } catch (const CalError &e) {
      fprintf(stderr, "spcat error: %s\n", e.what());
      return EXIT_FAILURE;
    }
    return 0;
  }

  // ── Legacy path (.var / .int) ─────────────────────────────────────────────

  std::unique_ptr<CalculationEngine> calc_engine =
      (engineType == "dpi")
      ? std::unique_ptr<CalculationEngine>(std::make_unique<DpiEngine>())
      : std::unique_ptr<CalculationEngine>(std::make_unique<SpinvEngine>());

  /* Open output stream */
  FILE *luout_f = file_helpers::open_output(fname[eout], "w");
  FileSink luout_sink(luout_f);

  /* Read input data from .int and .var files */
  CalCatInput input;
  if (!CalCatIO::readInput(fname[eint], fname[evar], input, calc_engine, &luout_sink)) {
    fclose(luout_f);
    return EXIT_FAILURE;
  }

  /* Open output streams based on iflg */
  FILE *luegy_f = luout_f;
  FILE *lustr_f = luout_f;
  int iflg = input.iflg;
  if (iflg >= 1000)
    iflg %= 1000;
  iflg %= 100;
  if (iflg >= 10) {
    lustr_f = file_helpers::open_output(fname[estr], "w");
    iflg %= 10;
  }
  if (iflg != 0) {
    luegy_f = file_helpers::open_output(fname[eegy], "w");
  }

  FileSink luegy_sink(luegy_f);
  FileSink lustr_sink(lustr_f);
  MemorySink lucat_sink; // buffered so we can sort before writing

  /* Create CalCat instance and run */
  CalCatOutput output;
  {
    SigintFlag sigint_guard; // (guard active for duration of run)
    CalCat calCat(calc_engine, &luout_sink, &lucat_sink, &luegy_sink, &lustr_sink);
    try
    {
      calCat.run(input, output);
    }
    catch (const CalError &e)
    {
      fprintf(stderr, "Catalog generation failed: %s\n", e.what());
      if (luegy_f != luout_f) fclose(luegy_f);
      if (lustr_f != luout_f) fclose(lustr_f);
      fclose(luout_f);
      return EXIT_FAILURE;
    }
  }

  /* Sort catalog lines then write to file */
  output.cat_lines = lucat_sink.drain_lines();
  output.sort_cat_lines();
  {
    FILE *fcat = file_helpers::open_output(fname[ecat], "w");
    for (const auto &line : output.cat_lines) {
      fputs(line.c_str(), fcat);
      fputc('\n', fcat);
    }
    fclose(fcat);
  }

  /* Close remaining output streams */
  if (luegy_f != luout_f)
    fclose(luegy_f);
  if (lustr_f != luout_f)
    fclose(lustr_f);
  fclose(luout_f);

  if (toml_out && !cat_toml.empty()) {
    try {
      save_cat_output_toml(output, cat_toml);
      fprintf(stderr, "spcat: wrote %s\n", cat_toml.c_str());
    } catch (const CalError &e) {
      fprintf(stderr, "spcat: --toml-out: failed to write %s: %s\n", cat_toml.c_str(), e.what());
    }
  }

  return 0;
}
