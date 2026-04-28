/***********************************************************************/
/*   Copyright (C) 1989, California Institute of Technology            */
/*   All rights reserved.  U. S. Government Sponsorship under          */
/*   NASA Contract NAS7-918 is acknowledged.                           */

/*   Herbert M. Pickett, 20 Mar 1989                                   */
/*   Revised version in c, 22 March 1999                               */
/*   Revised version in C++, 2026                                      */

/* THIS IS A GENERALIZED INTENSITY AND FREQUENCY CALCULATOR            */

/* THE HAMILTONIAN IS ASSUMED TO BE ORGANIZED INTO BLOCKS              */
/*     SUCH THAT ADJACENT SETS OF BLOCKS ARE CONNECTED BY TRANSITIONS  */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "engine/SpinvEngine.hpp"
#include "engine/DpiEngine.hpp"
#include "CalCat.hpp"
#include "CalCatIO.hpp"
#include "splib/calpgm_types.h"
#include "splib/ulib.h"
#include "spcat/OutputSink.hpp"
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
  std::unique_ptr<CalculationEngine> calc_engine = std::make_unique<SpinvEngine>();
  if (argc > 1) {
    if (strcasecmp(argv[1], "--dpi") == 0) {
      calc_engine = std::make_unique<DpiEngine>();
      char *name = *argv;
      argv++;
      *argv = name;
      argc--;
    } else if (strcasecmp(argv[1], "--spinv") == 0) {
      char *name = *argv;
      argv++;
      *argv = name;
      argc--;
    }
  }

  file_helpers::parse_file_args(argc, argv, NFILE, fname, ext);

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

  return 0;
}
