/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#ifndef CALFIT_HPP // Keep original include guard for CalFit.hpp
#define CALFIT_HPP

#include <cstdio>
#include <string>
#include <vector>
#include <memory>
#include "CalculationEngine.hpp" // Assuming bcd_t might be in a common header or here
// If bcd_t is from lsqfit.h or calpgm.h, ensure those are included before this file if needed by other .hpp files
// For now, let's assume bcd_t is known. If not, include lsqfit.h or calpgm.h here.
#include "lsqfit.h" // For bcd_t, MAXQN, etc.

#define NDCARD 130
#define LBLEN 10 // Length for parameter labels, as used in original calfit.c

/**
 * @brief Input data structure for CalFit class
 *
 * Holds the input parameters and line data for the fitting process.
 */
struct CalFitInput
{
  std::string title;
  std::vector<std::string> optionLines;

  int npar;   // Number of parameters
  int limlin; // Max requested lines from file, can be updated by linein
  int nitr;   // Number of iterations
  int nxpar_from_file; // The nxpar read from the second line of .par
  double marqp0;       // Marquardt parameter
  double xerrmx;         // Max (OBS-CALC)/ERROR
  double parfac_initial; // The parfac from the second line of .par
  double fqfacq;         // IR frequency scaling factor

  std::string namfil_from_options; // Name file for parameter labels
  int nfmt_from_options;  // Format for quantum numbers
  int itd_from_options;   // ITD option
  int ndbcd_from_options; // Number of BCD digits for parameters
  int noptn_count;

  int nfit;   // #/fitted parameters (return from getpar)
  int inpcor; // flag for input correlation/variance (return from getpar, updated by getvar)

  std::vector<unsigned char> idpar_data; // bcd_t is unsigned char
  std::vector<double> par_initial;
  std::vector<double> erp_initial;
  std::vector<char> parlbl_data_flat;
  std::vector<double> var_initial_from_getvar;

  std::vector<std::string> lineData_raw;
};

/**
 * @brief Output data structure for CalFit class
 *
 * Holds the results of the fitting process, including fitted parameters and statistics.
 */
struct CalFitOutput
{
  std::vector<double> par;   // Fitted parameters
  std::vector<double> erpar; // Estimated errors for parameters
  double xsqbest;            // Best RMS error after fitting
  int itr;                   // Number of iterations performed

  // Data needed for CalFitIO::writeOutput
  std::string title_for_output;
  std::vector<std::string> optionLines_for_output;
  int npar_final;
  int limlin_final;
  int nitr_final;             // actual iterations done
  int nxpar_final_for_header; // actual nxpar used for file output
  double marqlast_final;
  double xerrmx_final;
  double parfac_final_for_header; // parfac0 from original
  double fqfacq_final;
  int nfit_final;
  std::vector<unsigned char> idpar_final_for_output;
  std::vector<char> parlbl_final_for_output_flat;
  std::vector<double> erp_original_for_output; // erp from initial read for .par file
  std::vector<double> var_final_for_output;
  std::vector<double> dpar_final_for_putvar; // dpar array for putvar
};

/**
 * @brief Main class for spectroscopic parameter fitting
 *
 * Encapsulates the logic for fitting spectroscopic parameters to experimental line frequencies.
 * Uses a CalculationEngine (either SpinvEngine or DpiEngine) for underlying calculations.
 */
class CalFit
{
public:
  /**
   * @brief Constructor for CalFit
   * @param engineType Type of calculation engine to use ("spinv" or "dpi")
   * @param final_lufit_stream The file stream for the main .fit output log.
   */
  CalFit(const std::string &engineType, FILE *final_lufit_stream);

  /**
   * @brief Destructor for CalFit
   */
  ~CalFit();

  /**
   * @brief Run the fitting process
   * @param input Input data for fitting
   * @param output Output data for results (will be populated)
   * @return True if fitting is successful, false otherwise
   */
  bool run(const CalFitInput &input, CalFitOutput &output);

private:
  std::unique_ptr<CalculationEngine> calc;
  FILE *lufit; // Output file for fit results (owned by main, passed in)

  // Member variables for storing fitting data (C-style arrays)
  char *parlbl;
  double *par, *erp, *oldpar, *erpar, *dpar, *delbgn;
  double *fitbgn, *var, *oldfit, *fit;
  double *teig, *pmix;
  int *iperm;
  bcd_t *idpar; // bcd_t is unsigned char

  // Scalar state variables
  int m_npar; // Actual number of parameters after getpar
  int m_nfit;
  int ndfit, ndiag;
  long m_nsize_p;     // Changed from int to long based on maxmem return type if necessary (or keep size_t)
  int m_nxpar_actual; // Derived from input.nxpar_from_file and parameters
  int m_nxfit;
  int m_inpcor;
  double m_xerrmx;
  double m_parfac;  // Can change from input.parfac_initial (this is parfac0 in original main)
  double m_parfac0; // Store the initial parfac for output scaling logic
  double m_fqfacq;
  int m_ndbcd;
  int m_nfmt;
  int m_itd;
  int m_nline;  // Number of lines after linein
  int m_limlin; // Max lines requested
  int m_maxf_from_linein;
  int m_nblkpf_actual;
  int m_maxdm_actual;
  int m_ndfree0;
  double m_tiny = 1.5e-38; // from original main

  // Methods for breaking down the fitting process
  bool initializeParameters(const CalFitInput &input);
  bool processLinesAndSetupBlocks(const CalFitInput &input);
  bool performIteration(const CalFitInput &input, CalFitOutput &output);
  bool finalizeOutputData(const CalFitInput &input, CalFitOutput &output); // For populating CalFitOutput

  // Helper methods (kept as private, implementations copied from original if simple)
  int qnfmt2(int nqn, short *qnum, char *aqnum);
  int parer(double par_val, double errx, double dif, char *ptmp);  // Renamed par to par_val to avoid conflict
  int linein(FILE *luin, int *nline_val, int iqnfmt_val); // Renamed to avoid conflict if global linein exists
  int lineix(FILE *lu, int flg, int nline_val, int nblkpf_val, int iqnfmt_val);
  int getblk(int *iblk, int *indx, short *iqnum, int nblkpf_val, int ipos, int nqn_val);
};

#endif // CALFIT_HPP