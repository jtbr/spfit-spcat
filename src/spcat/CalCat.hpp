/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#ifndef CALCAT_HPP
#define CALCAT_HPP

#include <cstdio>
#include <string>
#include <vector>
#include <memory>
#include "engine/CalculationEngine.hpp"
#include "splib/lsqfit.h" // For bcd_t, MAXQN, BOOL, etc.
#include "common/Logger.hpp"
#include "spcat/OutputSink.hpp"

#define NTEMP 1001
#define TMAX 1000
#define MAXQNX 13

typedef struct {
  double *eigblk;
  double *egyblk;
  unsigned int nsizblk;
  int ixblk;
} SBLK;

/**
 * @brief Input data structure for CalCat class
 *
 * Holds all input data needed for catalog generation.
 */
struct CalCatInput
{
  // From .int file, line 1
  std::string title;

  // From .int file, line 2 (control parameters)
  int iflg;
  long itag;
  double qrot;
  int inblk;
  int lblk;
  double thrsh;
  double thrsh1;
  double fqmax;
  double tmq;
  int maxv;

  // Dipole data from .int file
  int ndip;
  std::vector<double> dip;
  std::vector<unsigned char> idip; // bcd_t
  std::vector<int> nvdip;
  std::vector<int> isimag;
  int npdip;

  // From .var file (via setopt + getpar + getvar)
  int npar;
  int nfit;
  int catqn;
  int nfmt;
  int itd;
  int ndbcd;
  std::vector<unsigned char> idpar;
  std::vector<double> par;
  std::vector<double> derv;
  std::vector<double> var;

  // From setfmt
  int nqn;
  std::vector<int> iqnfmtv;
};

/**
 * @brief Output data structure for CalCat class
 *
 * Holds summary results from catalog generation.
 * cat_lines/egy_lines/str_lines are populated when MemorySinks are drained
 * after run() (library/Python use); empty for FileSink (CLI) use.
 */
struct CalCatOutput
{
  long nline;
  double egymin;
  double qrot;    // initial Q value (for summary)
  BOOL ifdump;    // whether diagonalization was skipped

  int ntemp;
  double temp[NTEMP];
  double qsum[NTEMP];

  std::vector<std::string> cat_lines;
  std::vector<std::string> egy_lines;
  std::vector<std::string> str_lines;

  // Sort cat_lines by frequency (first 13 columns), matching the file-based
  // sortn step performed by the CLI after writing the .cat file.
  void sort_cat_lines();
};

/**
 * @brief Block cache state (replaces ibufof static variables)
 */
struct BlockCacheState
{
  FILE *scratch;
  long maxrec;
  long lsizb;
  int mempos;
  int orgpos;
  int nbsav;
  unsigned int maxdm;

  BlockCacheState() : scratch(nullptr), maxrec(0), lsizb(0),
                      mempos(0), orgpos(0), nbsav(0), maxdm(0) {}
};

/**
 * @brief Main class for spectroscopic catalog generation (SPCAT)
 *
 * Encapsulates the logic for computing spectroscopic line frequencies,
 * intensities, and partition functions from fitted Hamiltonian parameters.
 */
class CalCat
{
public:
  CalCat(std::unique_ptr<CalculationEngine> &calc_engine,
         OutputSink *luout, OutputSink *lucat, OutputSink *luegy, OutputSink *lustr,
         Logger &logger = Logger::defaultLogger());
  ~CalCat();

  void run(const CalCatInput &input, CalCatOutput &output);

private:
  std::unique_ptr<CalculationEngine> calc;
  Logger &m_logger;
  OutputSink *m_luout;
  OutputSink *m_lucat;
  OutputSink *m_luegy;
  OutputSink *m_lustr;

  // Phase methods (throw CalError on failure)
  void validateInput(const CalCatInput &input);
  void initializeParameters(const CalCatInput &input, CalCatOutput &output);
  void setupBlocks(const CalCatInput &input);
  void computeCatalog(const CalCatInput &input, CalCatOutput &output);
  void finalizeOutput(const CalCatInput &input, CalCatOutput &output);

  // Helper methods (in CalCat_helpers.cpp)
  static int qnfmt(short *iqu, int nqn, char *sqn);
  static int simt(int isiz, int jsiz, double *s, double *t, double *u, double *wk);
  SBLK *ibufof(int iblk, unsigned int ndel, SBLK *blk);
  static SBLK *sblk_alloc(int nstruct, unsigned mxdm);

  // Block cache state (replaces ibufof statics)
  BlockCacheState m_cache;

  // Control flags (derived from iflg)
  BOOL m_prir, m_prder, m_prfrq, m_preig, m_pregy, m_prstr;
  BOOL m_diag, m_ifdump;

  // Derived constants
  double m_fac, m_cmc, m_tmc, m_tmq, m_tmql;
  double m_fqmax, m_fqmin;
  double m_thrsh, m_thrsh1;
  double m_stcomp, m_strmn, m_starg;
  double m_bigerr, m_zero, m_telim;
  double m_qrot;

  // Block structure
  SBLK *m_blk;
  int m_nsav, m_nbsav;
  int m_nbkpj;
  unsigned int m_maxdm;
  unsigned int m_ndel;

  // Intensity matrices
  double **m_s;
  double *m_pmix;

  // Pointers to working copies of input data
  double *m_par;
  double *m_derv;
  double *m_var;
  bcd_t *m_idpar;
  double *m_dip;
  bcd_t *m_idip;
  int *m_nvdip;
  int *m_isimag;
  int *m_iqnfmtv;

  // Scalars from input
  int m_npar, m_nfit;
  int m_npdip, m_ndip;
  int m_nqn, m_nfmt;
  int m_itd;
  int m_catqn;
  long m_itag;
  int m_iflg;

  // QN format state
  int m_iqnfmt;
  int m_iposv;
  int m_maxv;
  int m_globfmt;
  int m_newfmt;
  int m_maxqn;

  // Block range
  int m_inblk, m_lblk;
};

#endif // CALCAT_HPP
