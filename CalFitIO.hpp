/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#ifndef CALFIT_IO_HPP
#define CALFIT_IO_HPP

#include <string>
#include <vector>
#include "CalFit.hpp" // Includes CalFitInput, CalFitOutput

/**
 * @brief Class for handling I/O operations for CalFit
 */
class CalFitIO
{
public:
  /**
   * @brief Static method to read input data from files
   * @param parFile Path to parameter file
   * @param linFile Path to line file
   * @param input Output parameter for input data (will be populated)
   * @param lufit_for_getpar File stream for getpar to write its log.
   * @return True if reading is successful, false otherwise
   */
  static bool readInput(const std::string &parFile, const std::string &linFile,
                        CalFitInput &input, FILE *lufit_for_getpar);

  /**
   * @brief Static method to write output data to files
   * @param par_filepath_final Path to the final .par output file (e.g. "base.par")
   * @param bak_filepath_original Path to the original .par file (now backup, e.g. "base.bak")
   * @param var_filepath_final Path to variance file (e.g. "base.var")
   * @param output Output data to write
   * @return True if writing is successful, false otherwise
   */
  static bool writeOutput(const std::string &par_filepath_final,
                          const std::string &bak_filepath_original,
                          const std::string &var_filepath_final,
                          const CalFitOutput &output,
                          const CalFitInput &original_input); // Added original_input for some fields
};

#endif // CALFIT_IO_HPP