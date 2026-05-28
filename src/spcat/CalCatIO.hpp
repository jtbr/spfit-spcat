/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#ifndef CALCAT_IO_HPP
#define CALCAT_IO_HPP

#include <string>
#include "CalCat.hpp"
#include "api/InputSchema.hpp"

/**
 * @brief Class for handling I/O operations for CalCat
 */
class CalCatIO
{
public:
  /**
   * @brief Read input data from .int and .var files
   * @param intFile Path to intensity/dipole file
   * @param varFile Path to variance/parameter file
   * @param input Output parameter for input data (will be fully populated)
   * @param calc_engine Calculation engine for setopt/setfmt calls
   * @param luout Output stream for diagnostic messages
   * @return True if reading is successful, false otherwise
   */
  static bool readInput(const std::string &intFile,
                        const std::string &varFile,
                        CalCatInput &input,
                        std::unique_ptr<CalculationEngine> &calc_engine,
                        OutputSink *luout);

  /**
   * @brief Write the catalog-settings preamble to a .out file (TOML path).
   *
   * Reproduces the header that CalCatIO::readInput writes when parsing legacy
   * .int / .var files: title line, ID/QSPINROT/QN line, log-strength line,
   * dipole list, .VAR FILE TITLE line, and PARAMETERS table.
   *
   * @param luout   Output file (mol.out, opened for writing)
   * @param ci      CatInput from load_cat_input_toml (ci.int_title = dipoles.toml title;
   *                ci.title = fitted.toml title; ci.dipoles = dipole list)
   * @param cci     CalCatInput from build_cat_input (has BCD-encoded idpar/idip, par, derv)
   */
  static void write_cat_preamble(FILE *luout, const CatInput &ci, const CalCatInput &cci);
};

#endif // CALCAT_IO_HPP
