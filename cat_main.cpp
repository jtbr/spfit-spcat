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
#include "SpinvEngine.hpp"
#include "DpiEngine.hpp"
#include "CalCat.hpp"
#include "CalCatIO.hpp"
#include "calpgm.h"

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

  filget(argc, argv, NFILE, fname, ext);

  /* Open output stream */
  FILE *luout = fopenq(fname[eout], "w");

  /* Read input data from .int and .var files */
  CalCatInput input;
  if (!CalCatIO::readInput(fname[eint], fname[evar], input, calc_engine, luout)) {
    fclose(luout);
    return EXIT_FAILURE;
  }

  /* Open output streams based on iflg */
  FILE *luegy = luout;
  FILE *lustr = luout;
  int iflg = input.iflg;
  if (iflg >= 1000)
    iflg %= 1000;
  iflg %= 100;
  if (iflg >= 10) {
    lustr = fopenq(fname[estr], "w");
    iflg %= 10;
  }
  if (iflg != 0) {
    luegy = fopenq(fname[eegy], "w");
  }

  /* Open catalog output stream */
  FILE *lucat = fopenq(fname[ecat], "w");

  /* Create CalCat instance and run */
  CalCatOutput output;
  {
    CalCat calCat(calc_engine, luout, lucat, luegy, lustr);
    calCat.run(input, output);
  }

  /* Close output streams */
  if (luegy != luout)
    fclose(luegy);
  if (lustr != luout)
    fclose(lustr);
  fclose(luout);
  fclose(lucat);

  /* Sort catalog file by frequency */
  sortn(fname[ecat], fname[ecat], FALSE);

  return 0;
}
