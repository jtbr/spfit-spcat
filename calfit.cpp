/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   25 March 1999: read option cards with fgetstr */
/*   30 Dec.  1999: include changes for dlsq */
/*   10 Oct.  2001: change fit diverging code */
/*   21 Sept. 2002: fix NRJ criterion */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint */

/**************************************************************************/
/*                                                                        */
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
//#include "calpgm.h"
#include "lsqfit.h"
#include "SpinvEngine.hpp"
#include "DpiEngine.hpp"
#include "CalFit.hpp"
#include "CalFitIO.hpp"

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
  static const char *ext[NFILE] = { "par", "lin", "fit", "bak", "var" };
  enum efile {epar, elin, efit, ebak, evar};
  char *fname[NFILE+1];
  std::string engineType = "spinv";

  // Choose the calculation engine based on the first argument
  if (argc > 1) {
    if (strcasecmp(argv[1], "--dpi") == 0) {
      engineType = "dpi";
      char *name = *argv; argv++; *argv = name; argc--; // Remove first argument
    } else if (strcasecmp(argv[1], "--spinv") == 0) {
      char *name = *argv; argv++; *argv = name; argc--; // Remove first argument
    }
  }

  // Open read and write files
  filget(argc, argv, NFILE, fname, ext);

  // Create CalFit instance with the selected engine
  CalFit calFit(engineType);

  // Prepare input and output structures
  CalFitInput input;
  CalFitOutput output;

  // Read input data using CalFitIO
  if (!CalFitIO::readInput(fname[epar], fname[elin], input)) {
    puts("Failed to read input files");
    return EXIT_FAILURE;
  }

  // Run the fitting process
  if (!calFit.run(input, output)) {
    puts("Fitting process failed");
    return EXIT_FAILURE;
  }

  // Write output data using CalFitIO
  if (!CalFitIO::writeOutput(fname[efit], fname[ebak], fname[evar], output)) {
    puts("Failed to write output files");
    return EXIT_FAILURE;
  }

  puts("FIT COMPLETE");
  return 0;
}                               /* MAIN */