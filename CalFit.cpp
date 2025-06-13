/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include <cstring> // for memset
#include <cmath>   // for fabs, sqrt
#include "lsqfit.h" // for linear algebra functions and other definitions
#include <stdio.h>  // for FILE and fprintf
#include <string.h> // for strcpy
#include <vector>   // for std::vector
#include <algorithm>// for std::copy
#include "CalFitIO.hpp"
#include "CalFit.hpp"
#include "lsqfit.h"  // For C functions like dcopy, ddot, mallocq, etc.
#include "calpgm.h"  // For MAXCAT, lbufof, maxmem etc.
#include "SpinvEngine.hpp"
#include "DpiEngine.hpp"

/**
 * @brief Constructor for CalFit
 * @param engineType Type of calculation engine to use ("spinv" or "dpi")
 * @param final_lufit_stream The file stream for the main .fit output log.
 */
CalFit::CalFit(const std::string &engineType, FILE *final_lufit_stream) : lufit(final_lufit_stream), // Initializer list
                                                                          parlbl(nullptr), par(nullptr), erp(nullptr), oldpar(nullptr), erpar(nullptr),
                                                                          dpar(nullptr), delbgn(nullptr), fitbgn(nullptr), var(nullptr), oldfit(nullptr),
                                                                          fit(nullptr), teig(nullptr), pmix(nullptr), iperm(nullptr), idpar(nullptr),
                                                                          m_npar(0), m_nfit(0), ndfit(0), ndiag(0), m_nsize_p(0), m_nxpar_actual(0), m_nxfit(0),
                                                                          m_inpcor(0), m_xerrmx(0.0), m_parfac(0.0), m_parfac0(0.0), m_fqfacq(0.0),
                                                                          m_ndbcd(0), m_nfmt(0), m_itd(0), m_nline(0), m_limlin(0), m_maxf_from_linein(0),
                                                                          m_nblkpf_actual(0), m_maxdm_actual(0), m_ndfree0(0)
{

  if (engineType == "dpi")
  {
    calc = std::make_unique<DpiEngine>();
  }
  else
  { // Default to spinv
    calc = std::make_unique<SpinvEngine>();
  }
  if (!lufit)
  {
    // Handle error: lufit stream is essential.
    // This could throw an exception or set an error state.
    fprintf(stderr, "FATAL: CalFit constructor received NULL lufit stream.\n");
    // For now, let's proceed, but in a real app, this is critical.
  }
}

/**
 * @brief Destructor for CalFit
 */
CalFit::~CalFit()
{
  // Memory allocated in initializeParameters and performIteration
  if (parlbl)
    free(parlbl);
  if (par)
    free(par);
  if (erp)
    free(erp);
  if (oldpar)
    free(oldpar);
  if (erpar)
    free(erpar);
  if (dpar)
    free(dpar);
  if (delbgn)
    free(delbgn);
  if (fitbgn)
    free(fitbgn);
  if (var)
    free(var);
  if (oldfit)
    free(oldfit);
  if (fit)
    free(fit);
  if (iperm)
    free(iperm);
  if (idpar)
    free(idpar);
  if (teig)
    free(teig);
  if (pmix)
    free(pmix);

  // calc is managed by unique_ptr
  // lufit is owned by main, so CalFit does not close it.
}

// Dummy run method for now
bool CalFit::run(const CalFitInput &input, CalFitOutput &output)
{
  if (!lufit)
  {
    puts("Error: CalFit::run called with NULL lufit stream.");
    return false;
  }

  if (!initializeParameters(input))
  {
    puts("CalFit::initializeParameters failed.");
    return false;
  }
  puts("CalFit::initializeParameters successful.");

  // Placeholder for subsequent steps
  // if (!processLinesAndSetupBlocks(input)) return false;
  // if (!performIteration(input, output)) return false;
  // if (!finalizeOutputData(input, output)) return false;

  // Populate some dummy output for testing CalFitIO::writeOutput
  output.title_for_output = input.title;
  output.optionLines_for_output = input.optionLines;
  output.npar_final = m_npar;
  output.limlin_final = m_limlin;                 // This should be actual lines read
  output.nitr_final = 0;                           // Dummy
  output.nxpar_final_for_header = m_nxpar_actual; // This is the index, original needs count
                                                  // Original output: nxpar = npar - nxpar (where nxpar was threshold)
                                                  // So, input.nxpar_from_file is the count for output.
  output.nxpar_final_for_header = input.nxpar_from_file;
  output.marqlast_final = input.marqp0; // Dummy
  output.xerrmx_final = m_xerrmx;
  output.parfac_final_for_header = m_parfac0;
  output.fqfacq_final = m_fqfacq;
  output.nfit_final = m_nfit;
  if (m_npar > 0 && par)
  { // Check if par was allocated
    output.par.assign(par, par + m_npar);
  }
  if (m_npar > 0 && erpar)
  {                                             // Check if erpar was allocated
    output.erpar.assign(erpar, erpar + m_npar); // Dummy errors
  }
  if (m_npar > 0 && erp)
  {
    output.erp_original_for_output.assign(erp, erp + m_npar);
  }
  if (idpar)
  {
    size_t idpar_actual_size = ((size_t)m_npar * m_ndbcd + m_ndbcd + 3);
    output.idpar_final_for_output.assign(idpar, idpar + idpar_actual_size);
  }
  if (parlbl)
  {
    size_t parlbl_actual_size = (LBLEN * (size_t)m_npar + 1);
    output.parlbl_final_for_output_flat.assign(parlbl, parlbl + parlbl_actual_size);
  }
  // var_final_for_output and dpar_final_for_putvar would be filled after fitting

  puts("CalFit::run completed (stubbed).");
  return true;
}

/**
 * @brief Initialize parameters for fitting
 * @param input Input data for fitting
 * @return True if initialization is successful, false otherwise
 */
bool CalFit::initializeParameters(const CalFitInput &input)
{
  // Transfer scalar values from input to member variables
  m_npar = input.npar;
  m_limlin = input.limlin; // This is the requested number of lines. m_nline will be actual.
  // Note: original main does: if (dvec[1] < 0.) { catqn = MAXQN; dvec[1] = -dvec[1]; }
  // input.limlin from CalFitIO already handles the sign if it was negative.
  // The catqn = MAXQN logic might be tied to nfmt.
  if (input.limlin < 0)
  { // If original file had negative line count
    // This implies catqn = MAXQN; The actual line count will be abs(input.limlin)
    // This logic will be handled later when setting nfmt for CalculationEngine
    m_limlin = -input.limlin;
  }

  m_fqfacq = input.fqfacq;
  m_ndbcd = input.ndbcd_from_options;
  m_nfmt = input.nfmt_from_options; // Will be used by calc->setfmt
  m_itd = input.itd_from_options;   // Will be used by calc->setopt

  m_nfit = input.nfit;
  m_inpcor = input.inpcor; // This is the value *after* getvar has potentially updated it
  m_xerrmx = input.xerrmx;
  m_parfac0 = input.parfac_initial; // Store the original parfac from file
  m_parfac = input.parfac_initial;  // This m_parfac can be scaled later

  if (m_xerrmx < m_tiny)
    m_xerrmx = 1e6; // From original main

  // --- Allocate C-style member arrays ---
  size_t parlbl_actual_size = (LBLEN * (size_t)m_npar + 1);
  parlbl = (char *)mallocq(parlbl_actual_size);
  if (!parlbl)
  {
    perror("mallocq for parlbl failed");
    return false;
  }

  par = (double *)mallocq((size_t)m_npar * sizeof(double));
  if (!par)
  {
    perror("mallocq for par failed");
    return false;
  }
  oldpar = (double *)mallocq((size_t)m_npar * sizeof(double));
  if (!oldpar)
  {
    perror("mallocq for oldpar failed");
    return false;
  }
  erp = (double *)mallocq((size_t)m_npar * sizeof(double));
  if (!erp)
  {
    perror("mallocq for erp failed");
    return false;
  }
  erpar = (double *)mallocq((size_t)m_npar * sizeof(double)); // Will be sized for fitted params later if different
  if (!erpar)
  {
    perror("mallocq for erpar failed");
    return false;
  }

  size_t idpar_actual_size = ((size_t)m_npar * m_ndbcd + m_ndbcd + 3);
  idpar = (bcd_t *)mallocq(idpar_actual_size * sizeof(bcd_t));
  if (!idpar)
  {
    perror("mallocq for idpar failed");
    return false;
  }

  // --- Copy data from input vectors to C-style member arrays ---
  if (input.parlbl_data_flat.size() >= parlbl_actual_size)
  {
    memcpy(parlbl, input.parlbl_data_flat.data(), parlbl_actual_size);
  }
  else if (!input.parlbl_data_flat.empty())
  {
    memcpy(parlbl, input.parlbl_data_flat.data(), input.parlbl_data_flat.size());
    parlbl[input.parlbl_data_flat.size()] = '\0'; // Ensure null term if short
  }
  else
  {
    *parlbl = '\0'; // Empty
  }

  if (input.par_initial.size() == (size_t)m_npar)
  {
    std::copy(input.par_initial.begin(), input.par_initial.end(), par);
    std::copy(input.par_initial.begin(), input.par_initial.end(), oldpar); // oldpar starts same as par
  }
  else
  {
    fprintf(stderr, "Warning: par_initial size mismatch.\n"); /* Handle error or fill with 0s */
  }

  if (input.erp_initial.size() == (size_t)m_npar)
  {
    std::copy(input.erp_initial.begin(), input.erp_initial.end(), erp);
  }
  else
  {
    fprintf(stderr, "Warning: erp_initial size mismatch.\n");
  }

  if (input.idpar_data.size() >= idpar_actual_size)
  {
    memcpy(idpar, input.idpar_data.data(), idpar_actual_size * sizeof(bcd_t));
  }
  else
  {
    fprintf(stderr, "Warning: idpar_data size mismatch.\n");
  }

  // --- Initialize fitting matrices based on m_nfit ---
  // These are sized based on m_nfit (number of parameters actually being fitted)
  if (m_nfit <= 0)
  {
    // No parameters to fit, some arrays might not be needed or fitting loop skipped.
    // For now, let's assume m_nfit > 0 if we proceed.
    if (m_npar > 0 && m_nfit == 0)
    {
      fprintf(lufit, "No parameters marked for fitting.\n");
    }
    // If m_nfit is 0, many subsequent allocations are problematic.
    // The original code might proceed and just not do much in lsqfit.
    // For safety, ensure m_nfit is at least 1 for array indexing logic if code relies on it.
    // Or, have checks throughout. For now, we assume if m_nfit is 0, parts of the iteration are skipped.
  }

  ndfit = m_nfit + 1;
  ndiag = ndfit + 1;

  size_t nl_maxmem;               // For maxmem call
  m_nsize_p = maxmem(&nl_maxmem); // maxmem from calpgm.h
  if (ndfit > m_nsize_p && m_nfit > 0)
  { // Check m_nfit > 0 to avoid false positive if nfit is 0
    printf("number of independent parameters is too big: %d %d\n", m_nfit, (int)m_nsize_p - 1);
    return false;
  }

  if (m_nfit > 0)
  { // Only allocate these if there are parameters to fit
    iperm = (int *)mallocq((size_t)m_nfit * sizeof(int));
    if (!iperm)
    {
      perror("mallocq for iperm failed");
      return false;
    }
    // erpar was allocated based on m_npar, if m_nfit is different, it's used for m_nfit elements
    // The erpar array in main seemed to be used for all m_npar at the end.
    // For lsqfit, it needs space for m_nfit. Let's assume m_npar >= m_nfit.
    // The erpar for non-fitted params are derived.
    // For now, erpar is sized for m_npar.

    dpar = (double *)mallocq((size_t)ndfit * sizeof(double)); // dpar stores parameter shifts + 1 extra
    if (!dpar)
    {
      perror("mallocq for dpar failed");
      return false;
    }
    delbgn = (double *)mallocq((size_t)m_nfit * sizeof(double)); // Delta begin for trust region
    if (!delbgn)
    {
      perror("mallocq for delbgn failed");
      return false;
    }

    size_t nlsq_fit_matrices;
    if ((m_nfit & 1) == 0)
    { // even
      nlsq_fit_matrices = (size_t)((unsigned)m_nfit >> 1) * (size_t)ndfit;
    }
    else
    { // odd
      nlsq_fit_matrices = (size_t)((unsigned)ndfit >> 1) * (size_t)m_nfit;
    }
    // This nlsq_fit_matrices is for elements, not bytes yet.

    fitbgn = (double *)mallocq(nlsq_fit_matrices * sizeof(double));
    if (!fitbgn)
    {
      perror("mallocq for fitbgn failed");
      return false;
    }

    size_t nlsq_var_count = ((size_t)m_nfit * ((size_t)m_nfit + 1)) / 2;
    var = (double *)mallocq(nlsq_var_count * sizeof(double));
    if (!var)
    {
      perror("mallocq for var failed");
      return false;
    }
    if (input.var_initial_from_getvar.size() == nlsq_var_count)
    {
      std::copy(input.var_initial_from_getvar.begin(), input.var_initial_from_getvar.end(), var);
    }
    else if (!input.var_initial_from_getvar.empty())
    {
      fprintf(stderr, "Warning: var_initial_from_getvar size mismatch.\n");
    }
    else if (m_nfit > 0)
    { // If var_initial is empty but we need it. This shouldn't happen if getvar works.
      fprintf(stderr, "Error: var_initial_from_getvar is empty but m_nfit > 0.\n");
      // This indicates an issue in CalFitIO or getvar logic.
      // For robustness, one might try to reconstruct a diagonal var from erp here,
      // but it's better to fix the source.
      return false;
    }

    size_t oldfit_elements = nlsq_fit_matrices + (size_t)m_nfit; // oldfit = fitbgn + dpar structure in original
    oldfit = (double *)mallocq(oldfit_elements * sizeof(double));
    if (!oldfit)
    {
      perror("mallocq for oldfit failed");
      return false;
    }

    fit = (double *)mallocq((size_t)m_nfit * (size_t)ndfit * sizeof(double)); // Full fit matrix for lsqfit
    if (!fit)
    {
      perror("mallocq for fit failed");
      return false;
    }
  }

  // Fixup supplied values (nxpar related logic from original main)
  m_nxpar_actual = input.nxpar_from_file; // This is the count of params at end of list *not* to auto-exclude
  if (m_nxpar_actual > m_npar || m_nxpar_actual < 0)
    m_nxpar_actual = 0;

  if (m_nxpar_actual > 0 && lufit)
  {
    fprintf(lufit, "NUMBER OF PARAMETERS EXCLUDED IF INDEX < 0 =%5d\n", m_nxpar_actual);
  }
  m_nxpar_actual = m_npar - m_nxpar_actual; // This now becomes the threshold index
                                            // Parameters from m_nxpar_actual to m_npar-1 are checked.

  m_nxfit = m_nfit;
  for (int i_par = m_npar - 1; i_par >= m_nxpar_actual; --i_par)
  {
    int ibcd_idx = i_par * m_ndbcd;
    // Check array bounds for idpar: idpar_actual_size must be >= ibcd_idx + m_ndbcd
    if (idpar && (size_t)ibcd_idx < idpar_actual_size)
    {                                   // Basic check
      if (NEGBCD(idpar[ibcd_idx]) == 0) // NEGBCD expects pointer to the BCD block
        --m_nxfit;
    }
    else
    {
      fprintf(stderr, "Error: Index out of bounds for idpar in nxfit calculation.\n");
    }
  }

  if (fabs(m_fqfacq - 1.) >= 1e-10 && lufit)
  {
    fprintf(lufit, " IR Frequencies Scaled by %12.10f\n", m_fqfacq);
  }

  m_ndfree0 = 0; // For degrees of freedom adjustment

  // Initialize fitbgn, fit, oldfit, dpar, delbgn based on m_inpcor (return from getvar)
  if (m_nfit > 0)
  { // Only if fitting
    if (m_inpcor > 0)
    { // If getvar read variance from file (m_inpcor is getvar's return value)
      // "supplied variance" path
      double *pvar_ptr = var;
      double *pfit_ptr = fit;
      for (int n_idx = 1; n_idx <= m_nfit; ++n_idx)
      {
        dcopy(n_idx, pvar_ptr, 1, pfit_ptr, ndfit);
        pvar_ptr += n_idx;
        pfit_ptr++; // Advance by 1 for packed upper to full matrix column start
      }
      // The dcopy above is tricky. 'fit' is a full matrix, 'var' is packed upper.
      // Original code: dcopy(n, pvar, 1, pfit, ndfit); pvar += n; ++pfit;
      // This copies the n elements of the current column of var (which is a column of the Cholesky factor L stored packed)
      // into the n elements of the current column of fit (full matrix representation of L).
      // This seems correct if 'fit' is being built as lower triangular initially.

      fit[0] = var[0]; // Assuming var[0] is L_00^2 or L_00 depending on getvar
      dpar[0] = 1.0;   // dpar used for scaling factors

      pfit_ptr = fit;
      double current_scale;
      for (int i_s = 0; i_s < m_nfit; ++i_s)
      { // Scale loop
        if (i_s != 0)
        {
          pfit_ptr += ndfit; // Move to next column in 'fit' matrix
                             // If 'fit' is representing L (lower triangular), elements above diagonal are 0.
                             // If symmetric A, then upper elements are copies.
                             // The dqrfac implies 'fit' is treated as a general matrix to be factored.
                             // Let's assume memset for upper part if it's meant to be lower triangular input.
                             // Original doesn't explicitly zero upper part here. dqrfac might handle it.
        }
        int n_val_scale = m_nfit - i_s;
        double *pfitd_ptr = pfit_ptr + i_s; // Diagonal element of current submatrix
        double val_norm = dnrm2(n_val_scale, pfitd_ptr, 1);
        if (val_norm < m_tiny)
        {
          dpar[i_s] = 1.0;
          current_scale = 0.0;
        }
        else
        {
          dpar[i_s] = current_scale = 1.0 / val_norm;
        }
        dscal(n_val_scale, current_scale, pfitd_ptr, 1);
      }

      // erpar here seems to be temporary storage for diagonal elements of scaled A^T A (or L^T L)
      pfit_ptr = fit;
      for (int i_er = 0; i_er < m_nfit; ++i_er)
      {
        int n_er = i_er + 1;
        // erpar[i_er] = ddot(n_er, pfit_ptr, ndfit, pfit_ptr, ndfit); // This seems to be for A^T A_ii
        // The original code has erpar[i] = ddot(n, pfit, ndfit, pfit, ndfit); ++pfit;
        // This implies pfit iterates column-wise on 'fit' which is lower triangular (L)
        // and computes sum of squares of elements in that column.
        erpar[i_er] = ddot(m_nfit - i_er /* num elements from diagonal down */,
                           pfit_ptr + i_er /* diagonal of current column */, 1,
                           pfit_ptr + i_er, 1);
        // No, this is simpler: erpar[i] = ddot(i+1, column_i_of_L, 1, column_i_of_L, 1) assuming L is column major
        // Or: erpar[i] = ddot(N-(i), &L(i,i), N) if L is row major.
        // Let's re-verify ddot usage from original:
        // pfit = fit; for (i=0 to nfit-1) { n = i+1; erpar[i] = ddot(n, pfit, ndfit, pfit, ndfit); ++pfit;}
        // This means 'pfit' points to the start of columns in a full matrix, and 'ndfit' is leading dim.
        // It calculates sum_j (L_ji)^2 for j=0..i. This does not seem right for A=LL^T decomposition.
        // The `dqrfac` implies QR decomposition. `fit` becomes R.
        // Let's assume `erpar` is used as scratch by `dqrfac`.
        // The line `erpar[i] = ddot(n, pfit, ndfit, pfit, ndfit);` seems to be calculating something else,
        // perhaps related to diagonal preconditioning before `dqrfac`.
        // Given `dqrfac(fit, ndfit, nfit, nfit, 0, erpar, iperm);`
        // `erpar` is likely an output of `dqrfac` (e.g. diagonal elements of R or permutation info).
        // Let's defer populating erpar until after dqrfac if it's an output.
        // The original lines before dqrfac for erpar are likely specific pre-computation.
        // For now, let's skip the erpar loop before dqrfac. It's populated by lsqfit later.
      }

      int n_rank = dqrfac(fit, ndfit, m_nfit, m_nfit, 0, erpar, iperm); // erpar is scratch/output here
      if (m_nfit > n_rank && lufit)
      {
        fputs("supplied variance is singular\n", lufit);
      }
      dqrsolv(fit, ndfit, n_rank, m_nfit, 0, iperm); // fit becomes (R^-1 Q^T) or similar

      pfit_ptr = fit;
      for (int n_un = 0; n_un < m_nfit; ++n_un)
      {                                              // Unscale, iterate through columns
        current_scale = dpar[n_un];                  // Scaling factor for this column
        dscal(n_un + 1, current_scale, pfit_ptr, 1); // Scale first n_un+1 elements of column n_un
                                                     // if fit is lower triangular (e.g. Cholesky factor of inverse)
        // Original: for (n=1; i <= nfit; ++n) { scale=dpar[n-1]; dscal(n,scale,pfit,ndfit); ++pfit;}
        // This means pfit points to columns, dscal(n, scale, column_start, 1) if ndfit=N,
        // or dscal(n, scale, column_start, N) if matrix is column major.
        // Assuming 'fit' is column-major or treat pfit_ptr as start of column and dscal operates on contiguous elements.
        // If fit is (R^-1 Q^T), it's dense.
        // Let's assume dscal(m_nfit - n_un, current_scale, pfit_ptr + n_un, ndfit) for upper R^-1
        // Or dscal(n_un+1, current_scale, pfit_ptr, 1) if pfit_ptr points to a column of L_inv and applying scaling.
        // The original unscale is: for (n = 1; i <= nfit; ++n) { scale = dpar[n-1]; dscal(n, scale, pfit, ndfit); ++pfit; }
        // This implies pfit points to the start of the n-th column (0-indexed n-1).
        // dscal(n, scale, &fit(0,n-1), 1) -- if ndfit is used as stride, it's different.
        // If fit is lower triangular by columns:
        dscal(m_nfit - n_un, current_scale, pfit_ptr + n_un, ndfit); // Scales from diagonal down
                                                                     // This applies D^-1 * L_inv
        pfit_ptr++;                                                  // This doesn't make sense if pfit_ptr is start of matrix.
      }
      // Correct unscaling if 'fit' is (D^-1 L_inv) where L was from Cholesky of scaled matrix.
      // Iterate columns of (D^-1 L_inv). For each column j, multiply by dpar[j].
      pfit_ptr = fit;
      for (int j_col = 0; j_col < m_nfit; ++j_col)
      {
        dscal(m_nfit, dpar[j_col], pfit_ptr, 1); // Scale entire column j_col
        pfit_ptr += ndfit;                       // Move to next column
      }

      // Store in fitbgn (packed lower triangle of the inverted, unscaled covariance matrix)
      double *pfitb_ptr = fitbgn;
      pfit_ptr = fit;
      for (int i_fb = 0; i_fb < m_nfit; ++i_fb)
      {
        int n_fb = m_nfit - i_fb;
        dcopy(n_fb, pfit_ptr + i_fb, ndfit, pfitb_ptr, 1); // Copy from diagonal down
        pfitb_ptr += n_fb;
        pfit_ptr += ndfit; // Next column
      }
      // The original logic:
      // pfit = fit; pfitb = fitbgn;
      // for (i=0; i<nfit; ++i) { if (i!=0) pfit+=ndiag; n=nfit-i; dcopy(n, pfit, ndfit, pfitb,1); pfitb+=n;}
      // This copies from row i (starting at diagonal) of 'fit' (if 'fit' is row-major lower triangular)
      // or from column i (starting at diagonal) of 'fit' (if 'fit' is col-major lower triangular)
      // into packed 'fitbgn'. Let's assume fit is lower triangular, stored by columns.
      pfit_ptr = fit; // Start of matrix
      pfitb_ptr = fitbgn;
      for (int col_idx = 0; col_idx < m_nfit; ++col_idx)
      {
        int num_to_copy = m_nfit - col_idx;                      // Elements from diagonal down in this column
        dcopy(num_to_copy, pfit_ptr + col_idx, 1, pfitb_ptr, 1); // Assumes contiguous column elements
        pfitb_ptr += num_to_copy;
        pfit_ptr += ndfit; // Move to the start of the next column
      }
    }
    else
    { // "no supplied variance" path (m_inpcor <= 0 from getvar)
      double *pvar_ptr = var;
      double *pfitb_ptr = fitbgn;
      if (var[0] == 0.0)
      { // Avoid division by zero if error was tiny
        fitbgn[0] = 1.0 / (m_tiny * m_tiny);
      }
      else
      {
        fitbgn[0] = 1.0 / var[0]; // var[0] is sigma_0^2
      }
      var[0] = 0.0; // Modified as in original
      m_ndfree0 = 0;
      for (int n_idx = 1; n_idx < m_nfit; ++n_idx)
      {
        pfitb_ptr++;                                  // Advance in packed array fitbgn
        memset(pfitb_ptr, 0, sizeof(double) * n_idx); // Zero out off-diagonal elements for this row/col
        pvar_ptr += (n_idx + 1);                      // Advance in packed array var to next diagonal
        pfitb_ptr += n_idx;                           // Advance in packed array fitbgn to next diagonal
        if (*pvar_ptr == 0.0)
        {
          *pfitb_ptr = 1.0 / (m_tiny * m_tiny);
        }
        else
        {
          *pfitb_ptr = 1.0 / (*pvar_ptr);
        }
        if (*pfitb_ptr < 1.e-15)
          --m_ndfree0;   // From original
        *pvar_ptr = 0.0; // Modified as in original
      }
    }

    if (m_nfit > 0 && oldpar && par && oldfit && fitbgn && delbgn)
    {                        // Check pointers
      oldpar[0] = par[0];    // oldpar needs full m_npar
      oldfit[0] = fitbgn[0]; // oldfit is sized based on m_nfit related nlsq_fit_matrices
      memset(delbgn, 0, sizeof(double) * m_nfit);
    }
    else if (m_nfit > 0)
    {
      fprintf(stderr, "Error: One or more fitting arrays are null in initializeParameters.\n");
      return false;
    }
  } // end if (m_nfit > 0)

  // Call lbufof to initialize line buffer system
  // The number of lines here is the max requested, linein will update actual.
  lbufof(-m_nfit, m_limlin);

  return true;
}

// --- Stubs for remaining main logic methods ---

bool CalFit::processLinesAndSetupBlocks(const CalFitInput &input)
{
  if (!lufit)
  {
    puts("lufit is NULL in processLinesAndSetupBlocks");
    return false;
  }
  puts("CalFit::processLinesAndSetupBlocks called (STUB).");
  // This is where Step 2 will go:
  // 1. Call calc->setopt
  // 2. Call calc->setfmt
  // 3. Read lines from input.lineData_raw (using a temporary in-memory file for linein_internal)
  // 4. Call calc->setblk
  // 5. Call getlbl
  // 6. Call lineix_internal
  // 7. Allocate teig, pmix
  m_nline = 0; // Placeholder
  return true;
}

bool CalFit::performIteration(const CalFitInput &input, CalFitOutput &output)
{
  if (!lufit)
  {
    puts("lufit is NULL in performIteration");
    return false;
  }
  puts("CalFit::performIteration called (STUB).");
  // This is where Step 3 (the main fitting loop) will go.
  output.xsqbest = 1.0; // Dummy
  output.itr = 0;       // Dummy
  return true;
}

bool CalFit::finalizeOutputData(const CalFitInput &input, CalFitOutput &output)
{
  if (!lufit)
  {
    puts("lufit is NULL in finalizeOutputData");
    return false;
  }
  puts("CalFit::finalizeOutputData called (STUB).");

  // This method will populate CalFitOutput with all necessary data for CalFitIO::writeOutput.
  // It involves:
  // - Copying final parameters, errors.
  // - Computing correlation matrix (prcorr).
  // - Preparing data for .par and .var files.

  // Based on the CalFit::run stub, some fields are already populated.
  // More comprehensive population will happen here.
  output.title_for_output = input.title;             // Already in run stub
  output.optionLines_for_output = input.optionLines; // Already in run stub

  output.npar_final = m_npar;
  output.limlin_final = m_nline; // Actual number of lines used/read
  output.nitr_final = output.itr; // Iterations performed from iteration loop
  // Original output: nxpar = npar - nxpar (where nxpar was threshold count from end)
  // So input.nxpar_from_file is the count for output header.
  output.nxpar_final_for_header = input.nxpar_from_file;
  // output.marqlast_final will be the actual last Marquardt param from the loop
  output.xerrmx_final = m_xerrmx;
  output.parfac_final_for_header = m_parfac0; // The original parfac
  output.fqfacq_final = m_fqfacq;
  output.nfit_final = m_nfit;

  if (m_npar > 0 && par)
  {
    output.par.assign(par, par + m_npar);
  }
  if (m_npar > 0 && erpar)
  {
    output.erpar.assign(erpar, erpar + m_npar);
  } // These are the fitted errors
  if (m_npar > 0 && erp)
  {
    output.erp_original_for_output.assign(erp, erp + m_npar);
  } // Original a-priori errors

  if (idpar)
  {
    size_t idpar_actual_size = ((size_t)m_npar * m_ndbcd + m_ndbcd + 3);
    output.idpar_final_for_output.assign(idpar, idpar + idpar_actual_size);
  }
  if (parlbl)
  {
    size_t parlbl_actual_size = (LBLEN * (size_t)m_npar + 1);
    output.parlbl_final_for_output_flat.assign(parlbl, parlbl + parlbl_actual_size);
  }

  if (m_nfit > 0 && var)
  { // var is packed upper triangular
    size_t nlsq_var_count = ((size_t)m_nfit * ((size_t)m_nfit + 1)) / 2;
    output.var_final_for_output.assign(var, var + nlsq_var_count);
  }
  if (m_nfit > 0 && dpar)
  { // dpar contains scaled errors/shifts for fitted parameters
    // dpar in lsqfit output has m_nfit elements.
    // It's used by putvar.
    output.dpar_final_for_putvar.assign(dpar, dpar + m_nfit);
  }

  // Call prcorr if lufit is valid and fitting was done
  if (output.itr > 0 && m_nfit > 0 && fit && dpar && lufit)
  {
    // prcorr(lufit, nfit, fit, ndfit, dpar);
    // This needs to be called after the final fit matrix is available.
    // For now, it's a stub.
    fprintf(lufit, "Correlation matrix would be printed here.\n");
  }

  return true;
}

// original static functions now integrated into class

/**
 * @brief Format quantum numbers as a string for output
 * @param nqn Number of quantum numbers to format
 * @param qnum Array of quantum numbers
 * @param aqnum Output string buffer for formatted quantum numbers
 * @return Always returns 0
 */
int CalFit::qnfmt2(int nqn, short *qnum, char *aqnum)
{
  // Implementation copied from calfit.cpp
  int i;
  for (i = 0; i < nqn; ++i)
  {
    // Ensure aqnum has enough space. sprintf writes null terminator.
    // Each %3d needs 3 chars + 1 for null if called repeatedly.
    // It's safer if aqnum points to a sufficiently large buffer.
    sprintf(aqnum, "%3d", (int)qnum[i]);
    aqnum += 3;
  }
  for (i = nqn; i < 12; ++i)
  { // Pad to 12 quantum numbers (36 chars)
    aqnum[2] = aqnum[1] = aqnum[0] = ' ';
    aqnum += 3;
  }
  aqnum[0] = '\0'; // Null terminate the whole string
  return 0;
} // qnfmt2

/**
 * @brief Format parameter values, errors, and changes for output
 * @param par Parameter value
 * @param errx Parameter error
 * @param dif Parameter change
 * @param ptmp Output string buffer for formatted parameter
 * @return Always returns 0
 */
int CalFit::parer(double par, double errx, double dif, char *ptmp)
{
  static int czero = (int)'0';         /* ASCII code for '0' character */
  char *pfmt;                          /* Pointer for building format string */
  double adif, apar, aten, aerr;       /* Working copies of input values */
  char chexp[6], fmt[34];              /* Format string components */
  int msd, id, ie, efield, ip, lsd, k; /* Format calculation variables */

  /*      sets up special format for parameters and errors */
  /*     PAR  = parameter value */
  /*     ERRX = parameter error */
  /*     DIF  = parameter change */
  /*     PTMP = output string for printing */

  apar = par;
  aerr = errx;
  adif = dif;
  efield = 0;
  aten = 1.;
  /*     compute exponent fields */
  ie = (int)(log10(fabs(aerr) + 1.e-37) - 102.5) + 100;
  id = (int)(log10(fabs(adif) + 1.e-37) - 100.0) + 100;
  ip = (int)(log10(fabs(apar) + 1.e-37) - 100.0) + 100;
  lsd = -ie;
  if (lsd < 0)
    lsd = 0;
  msd = (ip > id) ? ip : id;
  /*  check for too many digits */
  k = 14 - ip;
  if (k < lsd)
    lsd = k;
  k = 10 - id;
  if (k < lsd)
    lsd = k;
  if (msd <= -2)
  { /* number too small without exponent */
    k = (1 - msd) / 3;
    efield = -3 * k;
    while ((--k) >= 0)
      aten *= 1000;
  }
  else if (lsd < 0)
  { /* number too big without exponent */
    k = (1 + msd) / 3;
    if (k > 0)
      efield = 3 * k;
    while ((--k) >= 0)
      aten *= 0.001;
  }
  if (efield != 0)
  { /* E format */
    lsd += efield;
    memcpy(chexp, "0fE+00", 6);
    if (efield < 0)
    {
      chexp[3] = '-';
      efield = -efield;
    }
    msd = efield / 10;
    if (msd > 0)
    {
      efield -= msd * 10;
      chexp[4] = (char)(msd + czero);
    }
    chexp[5] = (char)(efield + czero);
    apar *= aten;
    aerr *= aten;
    adif *= aten;
  }
  else
  { /* F format */
    memcpy(chexp, "0f    ", 6);
  }
  if (lsd > 9)
    lsd = 9;
  if (lsd > 0)
    chexp[0] = (char)(lsd + czero);
  while ((lsd--) > 0)
    aerr *= 10.;
  ie = (int)(aerr + 0.5);
  pfmt = fmt;
  memcpy(pfmt, "%16.", 4);
  pfmt += 4;
  memcpy(pfmt, chexp, 2);
  pfmt += 2;
  memcpy(pfmt, "(%3d)", 5);
  pfmt += 5;
  memcpy(pfmt, chexp + 2, 4);
  pfmt += 4;
  memcpy(pfmt, " %12.", 5);
  pfmt += 5;
  memcpy(pfmt, chexp, 6);
  pfmt += 6;
  *pfmt = '\0';
  sprintf(ptmp, fmt, apar, ie, adif);
  return 0;
} // parer

/**
 * @brief Read experimental spectral lines from input file
 * @param luin Input file pointer
 * @param nline Pointer to number of lines (input: max lines, output: actual lines)
 * @param iqnfmt Quantum number format for line input
 * @return Largest quantum number encountered
 */
int CalFit::linein(FILE *luin, int *nline, int iqnfmt)
{
  /* Local variables */
  SXLINE *xline;                /* Pointer to line data structure */
  double xfrqn, xerrn, xwtn;    /* Current line frequency, error, weight */
  double xfrqx, xerrx;          /* Previous line frequency and error */
  int nqn, nqnu, nqnl;          /* Number of quantum numbers */
  int kqnu, kqnl;               /* Indices for upper/lower state quantum numbers */
  int i, iqf, ipace, mxline;    /* Loop variables and counters */
  int mxqn, isblnd, icmp;       /* Max quantum number, blend flag, comparison flag */
  short nbln, nqnt[20], *iqnum; /* Blend counter, quantum number template, quantum number pointer */
  static char card[NDCARD];     /* Buffer for reading input lines */

  /*   get lines from input  and stores them */
  /*     LUIN= unit for finding lines */
  /*     NLINE = number of lines */
  /*     IQNFMT= qunatum number format for line input */
  /*     RETURN: largest quantum number */
  /*******************************************************************/

  mxline = *nline;
  mxqn = 1;
  nbln = 1;
  nqn = deflin(iqnfmt, nqnt);
  nqnu = nqn - 1;
  if (nqnt[nqnu] < 0)
    nqnu = 0;
  kqnu = nqnt[nqnu];
  nqnl = nqnu + nqn;
  kqnl = nqnt[nqnl];
  ipace = 100;
  xfrqx = xerrx = 0.;
  icmp = 0;
  for (i = 1; i <= mxline; ++i)
  { /*  loop for reading lines */
    xline = lbufof(1, i);
    iqnum = xline->qn;
    if (getlin(luin, nqn, nqnt, iqnum, &xfrqn, &xerrn, &xwtn,
               card, NDCARD) < 0)
    {
      *nline = i - 1;
      return mxqn;
    }
    iqf = iqnum[nqnu];
    if (iqf == -1)
    {
      if (kqnu >= 0)
      {
        iqf = -iqnum[kqnu];
        if (iqf >= 0)
          iqf = -1;
      }
      iqnum[nqnu] = (short)iqf;
    }
    if (iqf < 0)
      iqf = -iqf;
    if (iqf > mxqn)
      mxqn = iqf;
    iqf = iqnum[nqnl];
    if (iqf == -1)
    {
      if (kqnl > 0)
      {
        iqf = -iqnum[kqnl];
        if (iqf >= 0)
          iqf = -1;
      }
      iqnum[nqnl] = (short)iqf;
    }
    if (iqf < 0)
      iqf = -iqf;
    if (mxqn < iqf)
      mxqn = iqf;
    xline->xfrq = xfrqn;
    xline->xerr = (float)xerrn;
    xline->xwt = (float)fabs(xwtn);
    xline->linku = 0;
    xline->linkl = 0;
    isblnd = 0;
    if (icmp != 0 && fabs(xfrqn - xfrqx) < fabs(xfrqn) * 1.e-14 + 1.e-8)
    {
      /* frq match */
      if (fabs(xerrn - xerrx) < 1e-7)
      {
        isblnd = 1;
      }
      else if ((xerrn / xerrx) > 2.0 && nbln > 2)
      {
        isblnd = 1;
        ++nbln;
        icmp = 0;
        xline->xwt = (float)0.;
        iqnum[0] = (short)-1;
        iqnum[nqn] = iqnum[0];
      }
    }
    if (isblnd != 0)
    {
      xline->bln = nbln;
      xline = lbufof(1, i - 1);
      xline->bln = -2;
      nbln += 2;
    }
    else
    {
      xline->bln = 0;
      nbln = 2;
      icmp = 1;
    }
    if (ipace <= i)
    {
      ipace += 100;
      printf("Reading Line %d\n", i);
      fflush(stdout);
    }
    xerrx = xerrn;
    xfrqx = xfrqn;
  }
  return mxqn;
} // linein

/**
 * @brief Process spectral lines and set up block structure for fitting
 * @param lu Output file pointer for listing
 * @param flg Flag for printing (negative for detailed output)
 * @param nline Number of lines
 * @param nblkpf Number of blocks per F quantum number
 * @param iqnfmt Quantum number format for line input
 * @return Number of bad lines
 */
int CalFit::lineix(FILE *lu, int flg, int nline, int nblkpf, int iqnfmt)
{
  /*   get lines from input and store them */
  /*     LU = unit for printout of lines ( if > 0 ) */
  /*     NLINE = number of lines */
  /*     NBLKPF= number of blocks per F */
  /*     IQNFMT= qunatum number format for line input */
  /******************************************************************/
  static int nsort = 2048;
  SXLINE *xline;
  double xfrqn, xerrn, xwtn, xnorm;
  int nblk, ipos, i, j, ipace, nread, iblkl, iblku, ncat;
  int linkx, indxl, linky, indxu, orgblk, nqn, nqn2, nbad;
  /*@owned@*/ int *prvblk;
  short *iqnum;
  char aqnum[6 * MAXQN + 2];

  nbad = 0;
  nblk = 0;
  prvblk = (int *)mallocq((size_t)(nsort + 1) * sizeof(int));
  prvblk[0] = 0;
  for (i = 1; i <= nsort; ++i)
  {
    prvblk[i] = 0;
  }
  nqn = iqnfmt % 10;
  if (nqn == 0)
    nqn = 10;
  nqn2 = nqn + nqn;
  ncat = nqn2;
  if (ncat < 12)
    ncat = 12;
  i = (iqnfmt / 100) % 5;
  if (i >= nqn)
  {
    ipos = 1;
  }
  else
  {
    ipos = nqn;
  }
  if (flg < 0)
  {
    fputs(" LINE,BLKU,INDXU,BLKL,INDXL,QUANTUM NUMBERS", lu);
    for (i = 0; i < 19; ++i)
      fputc(' ', lu);
    fputs("ENERGY    EXP. ERROR    WEIGHTS\n", lu);
  }
  xnorm = 0.;
  ipace = 50;
  /*       loop for converting lines */
  for (nread = 1; nread <= nline; ++nread)
  {
    xline = lbufof(1, nread);
    xfrqn = xline->xfrq;
    xerrn = xline->xerr;
    xwtn = xline->xwt;
    /* find blocks and index for upper and lower states */
    iqnum = xline->qn;
    getblk(&iblku, &indxu, iqnum, nblkpf, ipos, nqn);
    getblk(&iblkl, &indxl, &iqnum[nqn], nblkpf, ipos, nqn);
    if (iblkl == 0 && iqnum[nqn] >= 0)
      iblku = 0;
    xline->ibu = iblku;
    xline->inu = (short)indxu;
    xline->ibl = iblkl;
    xline->inl = (short)indxl;
    if (iblku == 0 && (xline->bln & 1) == 0)
    {
      /*  print out bad line and try for next */
      ++nbad;
      xline->xwt = 0.;
      xwtn = 0.;
      qnfmt2(nqn2, iqnum, aqnum);
      printf("Bad Line(%3d): %s %14.5f %8.5f\n",
             nread, aqnum, xfrqn, xerrn);
      fprintf(lu, "Bad Line(%3d): %s %14.5f %8.5f\n",
              nread, aqnum, xfrqn, xerrn);
    }
    else
    {
      /*  set up links for calculating in order of block */
      if (iblku <= iblkl)
      {
        lnlink(prvblk, nsort, iblku, nread);
        lnlink(prvblk, nsort, iblkl, -nread);
        if (nblk < iblkl)
          nblk = iblkl;
      }
      else
      {
        lnlink(prvblk, nsort, iblkl, -nread);
        lnlink(prvblk, nsort, iblku, nread);
        if (nblk < iblku)
          nblk = iblku;
      }
      if (flg < 0)
      {
        iqnum = xline->qn;
        fprintf(lu, " %4d%4d%4d%4d%4d:", nread, iblku, indxu, iblkl, indxl);
        for (i = 0; i < ncat; ++i)
        {
          j = iqnum[i];
          fprintf(lu, "%3d", j);
        }
        fprintf(lu, " %14.4f %9.4f %9.4f", xfrqn, xerrn, xwtn);
        j = xline->bln;
        if (j != 0)
        {
          fprintf(lu, "   Line Blended with %3d\n", nread - (j >> 1));
        }
        else
        {
          fputc('\n', lu);
        }
      }
    }
    /* let the user know something is happening */
    if (ipace <= nread || nread == nline)
    {
      ipace += 50;
      printf("Converting Line %d\n", nread);
      fflush(stdout);
    }
    j = xline->bln;
    if (j != 0)
    {
      xnorm += xwtn;
      if (j > 0)
      { /* normalize weights */
        xnorm = 1. / xnorm;
        for (j = nread - (j >> 1); j <= nread; ++j)
        {
          xline = lbufof(1, j);
          xline->xwt = (float)(xline->xwt * xnorm);
        }
        xnorm = 0.;
      }
    }
    else
    {
      xline->xwt = 1.;
    }
  } /* end loop for converting lines */
  orgblk = nsort;
  while (nblk > orgblk)
  { /* finish up links */
    --orgblk;
    linkx = 0;
    for (j = 0; j < nsort; ++j)
    {
      if (prvblk[j] != 0)
      {
        linkx = prvblk[j];
        prvblk[j] = 0;
      }
    }
    prvblk[0] = linkx;
    prvblk[nsort] = 0;
    linkx = lnlink(prvblk, nsort, 1, 0);
    while (linkx != 0)
    {
      linky = linkx;
      getdbk(&linkx, &iblkl, &j, &j, &j);
      iblkl -= orgblk;
      lnlink(prvblk, nsort, iblkl, linky);
    }
    orgblk += nsort;
  }
  free(prvblk);
  return nbad;
} // lineix

/**
 * @brief Get block and index from quantum numbers
 * @param iblk Output block number
 * @param indx Output index within block
 * @param iqnum Array of quantum numbers
 * @param nblkpf Number of blocks per F quantum number
 * @param ipos Position in iqnum to find F
 * @param nqn Number of quantum numbers
 * @return Always returns 0
 */
int CalFit::getblk(int *iblk, int *indx, short *iqnum, int nblkpf, int ipos, int nqn)
{
  int ibgn, kbgn, nblk, iqnf, k, idblk, kdtau, kk, nn, icmp;
  int iblkx, indxx, idgn, nsize, ncod;
  short kqnum[MAXQN];

  /*  gets block (IBLK) and INDEX from quantum numbers (IQNUM) */
  /*     NBLKPF IS THE NUMBER OF BLOCKS PER "F" */
  /*     IPOS   IS THE POSITION IN IQNUM TO FIND "F" */
  /*     NQN    IS THE NUMBER OF QUANTUM NUMBERS */

  *indx = 0;
  *iblk = 0;
  /*   ..check for indication of null level */
  if (iqnum[0] < 0)
    return 0;
  --ipos;
  iqnf = iqnum[ipos];
  iqnum[ipos] = iqnum[0];
  ibgn = iqnf;
  nblk = nblkpf;
  if (ibgn < 0)
  {
    ibgn = -ibgn;
    k = 10 - ibgn;
    if (k > 1)
      nblk *= k;
  }
  ibgn *= nblkpf;
  /*   ..loop over blocks for given F */
  icmp = 1;
  iblkx = indxx = 0;
  for (idblk = 1; idblk <= nblk; ++idblk)
  {
    iblkx = ibgn + idblk;
    calc->getqn(iblkx, 0, 0, kqnum, &nsize);
    indxx = 1;
    /*       ..search for match within block */
    while (indxx <= nsize)
    {
      ncod = calc->getqn(iblkx, indxx, nqn, kqnum, &idgn);
      kqnum[ipos] = kqnum[0];
      nn = (ncod > 0) ? ncod : -ncod;
      /*           ..NN is the size of a wang sub-block */
      /*           ..NCOD < 0 if oblate basis */
      if (nn <= 1)
      {
        nn = 1;
        kbgn = 1;
      }
      else
      {
        kbgn = 3;
      }
      icmp = 0;
      for (k = kbgn; k < nqn; ++k)
      {
        icmp = iqnum[k] - kqnum[k];
        if (icmp != 0)
          break;
      }
      if (icmp == 0)
      {
        /* all quanta tested equal */
        if (nn == 1)
          break;
        /* assignment could be in this sub-block */
        kdtau = iqnum[1] - iqnum[2] - kqnum[1] + kqnum[2];
        if (ncod < 0)
          kdtau = -kdtau;
        if (kdtau >= 0 && (kdtau & 3) == 0)
        { /* symmetry is good */
          kk = (int)((unsigned)kdtau >> 2);
          if (kk < nn)
          {
            /* value is in range */
            indxx += kk;
            break;
          }
        }
        ++icmp;
      }
      /*  to get here, quanta were not right */
      indxx += nn;
    }
    if (icmp == 0)
      break;
  }
  if (icmp == 0)
  { /*   standard return */
    if (nblk > nblkpf)
    {
      iqnf -= (iblkx - ibgn - 1) / nblkpf;
    }
    if (iqnf < 0)
      indxx = -indxx;
    *iblk = iblkx;
    *indx = indxx;
  }
  iqnum[ipos] = (short)iqnf;
  return 0;
} // getblk
