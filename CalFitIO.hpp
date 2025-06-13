/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#ifndef CALFIT_IO_HPP
#define CALFIT_IO_HPP

#include <string>
#include <vector>
#include "CalFit.hpp"

/**
 * @brief Class for handling I/O operations for CalFit
 */
class CalFitIO {
public:
    /**
     * @brief Static method to read input data from files
     * @param parFile Path to parameter file
     * @param linFile Path to line file
     * @param input Output parameter for input data
     * @return True if reading is successful, false otherwise
     */
    static bool readInput(const std::string& parFile, const std::string& linFile, CalFitInput& input);

    /**
     * @brief Static method to write output data to files
     * @param fitFile Path to fit output file
     * @param bakFile Path to backup file
     * @param varFile Path to variance file
     * @param output Output data to write
     * @return True if writing is successful, false otherwise
     */
    static bool writeOutput(const std::string& fitFile, const std::string& bakFile, const std::string& varFile, const CalFitOutput& output);
};

#endif // CALFIT_IO_HPP
