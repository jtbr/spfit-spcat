/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#ifndef CALCAT_IO_HPP
#define CALCAT_IO_HPP

#include <string>
#include "CalCat.hpp"

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
};

#endif // CALCAT_IO_HPP
