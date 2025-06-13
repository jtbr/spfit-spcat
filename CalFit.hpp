/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#ifndef CALFIT_HPP
#define CALFIT_HPP

#include <string>
#include <vector>
#include <memory>
#include <cstdio>
#include "CalculationEngine.hpp"

#define MAXQN 10
#define NDCARD 130
#define LBLEN 16 // Length for parameter labels, as used in original calfit.c

/**
 * @brief Input data structure for CalFit class
 *
 * Holds the input parameters and line data for the fitting process.
 */
struct CalFitInput {
    int npar;              // Number of parameters
    int limlin;            // Number of lines requested
    int nitr;              // Number of iterations
    int nxpar;             // Number of parameters excluded if index < 0
    double marqp0;         // Marquardt parameter
    double xerrmx;         // Max (OBS-CALC)/ERROR
    double parfac;         // Parameter error scaling factor
    double fqfacq;         // IR frequency scaling factor
    std::string namfil;    // Name file for parameter labels
    int noptn;             // Number of option lines
    int nfmt;              // Format for quantum numbers
    int itd;               // ITD option
    int ndbcd;             // Number of BCD digits for parameters
    std::vector<std::string> optionLines; // Option lines from input file
    std::string title;     // Title from parameter file
    int nfit;              // Number of fitted parameters
    int inpcor;            // Flag for input correlation/variance
    // Additional input data for parameters and lines will be added as needed
};

/**
 * @brief Output data structure for CalFit class
 *
 * Holds the results of the fitting process, including fitted parameters and statistics.
 */
struct CalFitOutput {
    std::vector<double> par;    // Fitted parameters
    std::vector<double> erpar;  // Estimated errors for parameters
    double xsqbest;             // Best RMS error after fitting
    int itr;                    // Number of iterations performed
    // Additional output data will be added as needed
};

/**
 * @brief Utility class for handling I/O operations for CalFit
 */
class CalFitIO {
  public:
    /**
     * @brief Read input data from files
     * @param parFile Path to parameter file
     * @param linFile Path to line file
     * @param input Output parameter for input data
     * @return True if reading is successful, false otherwise
     */
    static bool readInput(const std::string& parFile, const std::string& linFile, CalFitInput& input);

    /**
     * @brief Write output data to files
     * @param fitFile Path to fit output file
     * @param bakFile Path to backup file
     * @param varFile Path to variance file
     * @param output Output data to write
     * @return True if writing is successful, false otherwise
     */
    static bool writeOutput(const std::string& fitFile, const std::string& bakFile, const std::string& varFile, const CalFitOutput& output);
};

/**
 * @brief Main class for spectroscopic parameter fitting
 *
 * Encapsulates the logic for fitting spectroscopic parameters to experimental line frequencies.
 * Uses a CalculationEngine (either SpinvEngine or DpiEngine) for underlying calculations.
 */
class CalFit {
  public:
    /**
     * @brief Constructor for CalFit
     * @param engineType Type of calculation engine to use ("spinv" or "dpi")
     */
    CalFit(const std::string& engineType = "spinv");

    /**
     * @brief Destructor for CalFit
     */
    ~CalFit();

    /**
     * @brief Run the fitting process
     * @param input Input data for fitting
     * @param output Output data for results
     * @return True if fitting is successful, false otherwise
     */
    bool run(const CalFitInput& input, CalFitOutput& output);

  private:
    std::unique_ptr<CalculationEngine> calc;

    // Member variables for storing fitting data
    char *parlbl;
    double *par, *erp, *oldpar, *erpar, *dpar, *delbgn;
    double *fitbgn, *var, *oldfit, *fit;
    int *iperm;
    bcd_t *idpar;
    FILE *lufit; // Output file for fit results
    int nfit, ndfit, ndiag, nsize_p;
    int nxpar, nxfit;
    int inpcor;
    double xerrmx, scale, tiny = 1.5e-38;
    int ndfree0;

    // Methods for breaking down the fitting process
    bool initializeParameters(const CalFitInput& input);
    bool processLines(const CalFitInput& input);
    bool performIteration(const CalFitInput& input, CalFitOutput& output);
    bool calculateErrors(const CalFitInput& input, CalFitOutput& output);

    // Helper methods
    int qnfmt2(int nqn, short *qnum, char *aqnum);
    int parer(double par, double errx, double dif, char *ptmp);
    int linein(FILE *luin, int *nline, int iqnfmt);
    int lineix(FILE *lu, int flg, int nline, int nblkpf, int iqnfmt);
    int getblk(int *iblk, int *indx, short *iqnum, int nblkpf, int ipos, int nqn);
};

#endif // CALFIT_HPP
