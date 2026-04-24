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
#include "Logger.hpp"

#define NDCARD 130
#define LBLEN 10 // Length for parameter labels, as used in original calfit.c

/**
 * @brief Input data structure for CalFit class
 *
 * Holds the input parameters and line data for the fitting process.
 */
struct CalFitInput
{
  // From .par file header/early lines
  std::string title;   // 1st line

  // 2nd line
  int npar;   // Number of parameters
  size_t limlin; // Max requested lines from file, can be updated by linein
  int nitr;            // Number of iterations
  int nxpar_from_file; // The nxpar read from the second line of .par
  double marqp0;       // Marquardt parameter
  double xerrmx;         // Max (OBS-CALC)/ERROR
  double parfac_initial; // The parfac from the second line of .par
  double fqfacq;         // IR frequency scaling factor
  int catqn;             // Quantum numbers to write

  std::vector<std::string> raw_option_lines_from_par;

  // From CalculationEngine::setopt (called by CalFitIO)
  int nfmt_cat_from_setopt; // Catalog-related nfmt output by setopt
  int itd_from_setopt;            // ICD option
  int ndbcd_from_setopt;          // Number of BCD digits for parameters
  std::string namfil_from_setopt; // Name file for parameter labels
  int noptn_read_by_setopt; // Number of option lines successfully processed by setopt

  // From getpar (called by CalFitIO)
  int nfit;   // #/fitted parameters (return from getpar)
  int inpcor; // flag for input correlation/variance (return from getpar, updated by getvar)
  std::vector<unsigned char> idpar_data; // bcd_t is unsigned char
  std::vector<double> par_initial;
  std::vector<double> erp_initial;
  std::vector<char> parlbl_data_flat;

  // From getvar (called by CalFitIO)
  std::vector<double> var_initial_from_getvar;

  // From .lin file
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
  int nitr_final_actual;  // actual iterations done
  int nxpar_for_header;   // actual nxpar used for file output
  double marqlast_final;
  double xerrmx_final;
  double parfac_for_header; // parfac0 from original
  double fqfacq_final;
  int nfit_final;
  std::vector<unsigned char> idpar_final_for_output; // this->idpar (m_npar * m_ndbcd + ...)
  std::vector<char> parlbl_final_for_output_flat; // this->parlbl
  std::vector<double> erp_original_for_output; // this->erp from initial read for .par file (a priori errors)
  std::vector<double> var_final_for_output;  // this->var (packed, m_nfit elements)
  std::vector<double> dpar_final_for_putvar; // this->dpar array for putvar ((m_nfit elements, from lsqfit's enorm, scaled for error)
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
   * @param calc_engine Calculation Engine implementation object (currently either SpinvEngine or DpiEngine)
   * @param final_lufit_stream The file stream for the main .fit output log.
   */
  CalFit(std::unique_ptr<CalculationEngine> &calc_engine, FILE *final_lufit_stream,
         Logger &logger = Logger::defaultLogger());

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
  void run(const CalFitInput &input, CalFitOutput &output);

private:
  std::unique_ptr<CalculationEngine> calc;
  FILE *lufit; // Output file for fit results (owned by main, passed in)
  Logger &m_logger;

  // Member variables for storing fitting data (C-style arrays)
  char *parlbl;
  double *par, *erp, *oldpar, *erpar, *dpar, *delbgn;
  double *fitbgn, *var, *oldfit, *fit;
  double *teig, *pmix;
  int *iperm;
  bcd_t *idpar; // bcd_t is unsigned char

  // Scalar state variables, set from CalFitInput
  int m_npar; // Actual number of parameters after getpar
  int m_nfit;
  int ndfit; // leading dimension of fit matrix = nfit+1; element (r,c) = fit[r + c*ndfit]
  int m_catqn;
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
  size_t m_limlin; // Max lines requested
  int m_nitr_requested; // #/ iterations requested
  int m_maxf_from_linein;
  int m_nblkpf_actual;
  int m_maxdm_actual;
  int m_ndfree0;
  double m_tiny = 1.5e-38; // from original main

  char m_namfil_buffer[NDCARD]; // To store namfil if needed, though input has it

  // State for iteration loop
  int m_itr; // Current iteration count
  double m_xsqbest;
  double m_marqlast; // Last Marquardt parameter for output
  double m_marqp[3]; // Marquardt parameters for lsqfit
//  double m_parfac0;  // Original parfac for scaling logic
  // m_parfac is already a member, potentially scaled
  int m_nqn_for_iteration;

  // Methods for breaking down the fitting process (throw CalError on failure)
  void validateInput(const CalFitInput &input);
  void initializeParameters(const CalFitInput &input);
  void processLinesAndSetupBlocks(const CalFitInput &input);
  void performIteration(const CalFitInput &input, CalFitOutput &output);
  void finalizeOutputData(const CalFitInput &input, CalFitOutput &output);

  // Helper methods (kept as private, implementations copied from original, now in CalFit_helpers.cpp)
  int qnfmt2(int nqn, short *qnum, char *aqnum);
  int parer(double par_val, double errx, double dif, char *ptmp);
  int linein(const std::vector<std::string> &lines, int *nline_val, int iqnfmt);
  int lineix(FILE *lu, int flg, int nline_val, int nblkpf_val, int iqnfmt_val);
  int getblk(int *iblk, int *indx, short *iqnum, int nblkpf_val, int ipos, int nqn_val);
};

#endif // CALFIT_HPP