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
  static const char *ext[NFILE] = {"par", "lin", "fit", "bak", "var"};
  enum efile
  {
    epar,
    elin,
    efit,
    ebak,
    evar
  };
  char *fname[NFILE + 1];
  std::string engineType = "spinv"; // Default

  // Choose the calculation engine based on the first argument
  if (argc > 1)
  {
    if (strcasecmp(argv[1], "--dpi") == 0)
    {
      engineType = "dpi";
      if (argc > 2)
      { // Ensure there's a base filename after --dpi
        char *name = *argv;
        argv++;
        *argv = name;
        argc--; // Remove --dpi
      }
      else
      {
        // Handle case where only --dpi is provided without filename.
        // For now, assume filget handles it or prints usage.
      }
    }
    else if (strcasecmp(argv[1], "--spinv") == 0)
    {
      // engineType is already spinv (default)
      if (argc > 2)
      {
        char *name = *argv;
        argv++;
        *argv = name;
        argc--; // Remove --spinv
      }
    }
  }

  std::unique_ptr<CalculationEngine> calc_engine;
  if (engineType == "dpi")
  {
    calc_engine = std::make_unique<DpiEngine>();
  }
  else
  {
    calc_engine = std::make_unique<SpinvEngine>();
  }
  printf("Using %s engine.\n", engineType.c_str());

  // Get filenames
  filget(argc, argv, NFILE, fname, ext); // fname[epar]="base.par", etc.

  // Backup .par file to .bak file
  // Original: filbak(fname[epar], fname[ebak]);
  // fname[ebak] is "base.bak"
  // This means "base.par" is copied to "base.bak"
  // The program then reads from "base.par" (which is now fname[ebak] in original logic if lubak was fname[ebak])
  // And writes the new .par to fname[epar].
  // Let's stick to: read from original fname[epar], backup fname[epar] to fname[ebak].
  // CalFitIO::writeOutput will write to a new fname[epar].

  if (!filbak(fname[epar], fname[ebak]))
  { // filbak returns 0 on success in some implementations TODO
    printf("Warning: Failed to create backup file %s from %s. Proceeding without backup.\n", fname[ebak], fname[epar]);
    // Or exit if backup is critical
  }
  else
  {
    printf("Backup of %s to %s created successfully.\n", fname[epar], fname[ebak]);
  }

  // Open the main .fit output stream
  FILE *lufit_stream = fopen(fname[efit], "w");
  if (!lufit_stream)
  {
    printf("Error: Unable to open fit output file '%s'.\n", fname[efit]);
    return EXIT_FAILURE;
  }
//  printf("Successfully opened fit output file '%s' for writing.\n", fname[efit]);

  // Create CalFit instance with the selected engine and lufit stream
  CalFit calFit(calc_engine, lufit_stream);

  // Prepare input and output structures
  CalFitInput input;
  CalFitOutput output;

  // Read input data using CalFitIO
  // We read from the original .par file (fname[epar])
  if (!CalFitIO::readInput(fname[ebak], fname[elin], input, calc_engine, lufit_stream))
  {
    printf("Failed to read input files using CalFitIO::readInput.\n");
    fclose(lufit_stream);
    return EXIT_FAILURE;
  }
  printf("CalFitIO::readInput completed.\n");

  // Run the fitting process
  if (!calFit.run(input, output))
  {
    printf("Fitting process failed in CalFit::run.\n");
    fclose(lufit_stream); // Close lufit before exiting on error
    return EXIT_FAILURE;
  }
  printf("CalFit::run completed.\n");

  // Close the main .fit output stream (owned by main)
  fclose(lufit_stream);
  printf("Closed fit output file '%s'.\n", fname[efit]);

  // Write output data using CalFitIO
  // CalFitIO::writeOutput writes to:
  // - fname[epar] (new parameter file)
  // - fname[evar] (new variance file)
  // fname[ebak] (the backup of the original .par) is kept.
  if (!CalFitIO::writeOutput(fname[epar], fname[ebak], fname[evar], output, input))
  {
    printf("Failed to write output files using CalFitIO::writeOutput.\n");
    return EXIT_FAILURE;
  }
  printf("CalFitIO::writeOutput completed.\n");

  puts("FIT COMPLETE (from main)");
  return 0;
} // main