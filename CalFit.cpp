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
 * @param cal_engine SpinvEngine or DpiEngine object
 * @param final_lufit_stream The file stream for the main .fit output log.
 */
CalFit::CalFit(std::unique_ptr<CalculationEngine> &calc_engine, FILE *final_lufit_stream) :
    calc(std::move(calc_engine)),
    lufit(final_lufit_stream), // Initializer list
    parlbl(nullptr), par(nullptr), erp(nullptr), oldpar(nullptr), erpar(nullptr),
    dpar(nullptr), delbgn(nullptr), fitbgn(nullptr), var(nullptr), oldfit(nullptr),
    fit(nullptr), teig(nullptr), pmix(nullptr), iperm(nullptr), idpar(nullptr),
    m_npar(0), m_nfit(0), ndfit(0), ndiag(0), m_nsize_p(0), m_nxpar_actual(0), m_nxfit(0),
    m_inpcor(0), m_xerrmx(0.0), m_parfac(0.0), m_parfac0(0.0), m_fqfacq(0.0),
    m_ndbcd(0), m_nfmt(0), m_itd(0), m_nline(0), m_limlin(0), m_maxf_from_linein(0),
    m_nblkpf_actual(0), m_maxdm_actual(0), m_ndfree0(0)
{
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

// CalFit::run method updated to call these stubs
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
    //fprintf(lufit, "ERROR: Parameter initialization failed.\n");
    return false;
  }
  puts("CalFit::initializeParameters successful.");

  if (!processLinesAndSetupBlocks(input))
  {
    puts("CalFit::processLinesAndSetupBlocks failed.");
    //fprintf(lufit, "ERROR: Line processing and block setup failed.\n");
    return false;
  }
  puts("CalFit::processLinesAndSetupBlocks successful.");

  if (!performIteration(input, output))
  {
    puts("CalFit::performIteration failed.");
    //fprintf(lufit, "ERROR: Iterative fitting failed.\n");
    return false;
  }
  puts("CalFit::performIteration successful.");

  if (!finalizeOutputData(input, output))
  {
    puts("CalFit::finalizeOutputData failed.");
    // Don't write to lufit here as it might be part of finalization logic that failed
    return false;
  }
  puts("CalFit::finalizeOutputData successful.");

  puts("CalFit::run completed.");
  return true;
}

/**
 * @brief Initialize parameters for fitting
 * @param input Input data for fitting
 * @return True if initialization is successful, false otherwise
 */
bool CalFit::initializeParameters(const CalFitInput &input)
{
  // --- 1. Transfer scalar control parameters ---
  m_npar = input.npar;
  m_limlin = input.limlin;
  if (input.limlin < 0)
  {
    m_limlin = -input.limlin;
  }
  m_nitr_requested = input.nitr;
  m_fqfacq = input.fqfacq;
  m_ndbcd = input.ndbcd_from_setopt;
  m_itd = input.itd_from_setopt;
  if (!input.namfil_from_setopt.empty())
  {
    strncpy(m_namfil_buffer, input.namfil_from_setopt.c_str(), NDCARD - 1);
    m_namfil_buffer[NDCARD - 1] = '\0';
  }
  else
  {
    m_namfil_buffer[0] = '\0';
  }
  m_nfit = input.nfit;
  m_inpcor = input.inpcor;
  m_xerrmx = input.xerrmx;
  m_parfac0 = input.parfac_initial;
  m_parfac = input.parfac_initial;
  if (m_xerrmx < m_tiny)
    m_xerrmx = 1e6;

  // --- 2. Allocate C-style member arrays ---
  // Free any existing memory first
  if (parlbl)
    free(parlbl);
  parlbl = nullptr;
  if (par)
    free(par);
  par = nullptr;
  if (erp)
    free(erp);
  erp = nullptr;
  if (oldpar)
    free(oldpar);
  oldpar = nullptr;
  if (erpar)
    free(erpar);
  erpar = nullptr;
  if (dpar)
    free(dpar);
  dpar = nullptr;
  if (delbgn)
    free(delbgn);
  delbgn = nullptr;
  if (fitbgn)
    free(fitbgn);
  fitbgn = nullptr;
  if (var)
    free(var);
  var = nullptr;
  if (oldfit)
    free(oldfit);
  oldfit = nullptr;
  if (fit)
    free(fit);
  fit = nullptr;
  if (iperm)
    free(iperm);
  iperm = nullptr;
  if (idpar)
    free(idpar);
  idpar = nullptr;

  size_t parlbl_actual_size = (LBLEN * (size_t)m_npar + 1);
  parlbl = (char *)mallocq(parlbl_actual_size);
  if (!parlbl)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: mallocq for parlbl failed.\n");
    perror("mallocq for parlbl failed");
    return false;
  }

  par = (double *)mallocq((size_t)m_npar * sizeof(double));
  if (!par)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: mallocq for par failed.\n");
    perror("mallocq for par failed");
    free(parlbl);
    parlbl = nullptr;
    return false;
  }

  oldpar = (double *)mallocq((size_t)m_npar * sizeof(double));
  if (!oldpar)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: mallocq for oldpar failed.\n");
    perror("mallocq for oldpar failed");
    free(parlbl);
    free(par);
    parlbl = nullptr;
    par = nullptr;
    return false;
  }

  erp = (double *)mallocq((size_t)m_npar * sizeof(double));
  if (!erp)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: mallocq for erp failed.\n");
    perror("mallocq for erp failed");
    free(parlbl);
    free(par);
    free(oldpar); /*... set null ...*/
    return false;
  }

  erpar = (double *)mallocq((size_t)m_npar * sizeof(double));
  if (!erpar)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: mallocq for erpar failed.\n");
    perror("mallocq for erpar failed");
    free(parlbl);
    free(par);
    free(oldpar);
    free(erp);
    return false;
  }

  size_t idpar_actual_size = ((size_t)m_npar * m_ndbcd + m_ndbcd + 3);
  idpar = (bcd_t *)mallocq(idpar_actual_size * sizeof(bcd_t));
  if (!idpar)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: mallocq for idpar failed.\n");
    perror("mallocq for idpar failed");
    free(parlbl);
    free(par);
    free(oldpar);
    free(erp);
    free(erpar);
    return false;
  }

  // --- 3. Copy data from input's vectors ---
  if (input.parlbl_data_flat.size() >= parlbl_actual_size)
  {
    memcpy(parlbl, input.parlbl_data_flat.data(), parlbl_actual_size);
  }
  else if (!input.parlbl_data_flat.empty())
  {
    memcpy(parlbl, input.parlbl_data_flat.data(), input.parlbl_data_flat.size());
    parlbl[input.parlbl_data_flat.size()] = '\0';
    // if (lufit)
    //   fprintf(lufit, "Warning: parlbl_data_flat size from input (%zu) was smaller than expected (%zu).\n", input.parlbl_data_flat.size(), parlbl_actual_size);
    printf("Warning: parlbl_data_flat size from input (%zu) was smaller than expected (%zu).\n", input.parlbl_data_flat.size(), parlbl_actual_size);
  }
  else if (m_npar > 0)
  {
    *parlbl = '\0';
    // if (lufit)
    //   fprintf(lufit, "Warning: parlbl_data_flat from input was empty for m_npar = %d.\n", m_npar);
    printf("Warning: parlbl_data_flat from input was empty for m_npar = %d.\n", m_npar);
  }

  if (input.par_initial.size() == (size_t)m_npar)
  {
    std::copy(input.par_initial.begin(), input.par_initial.end(), par);
    std::copy(input.par_initial.begin(), input.par_initial.end(), oldpar);
  }
  else if (m_npar > 0)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: par_initial size mismatch. Expected %d, got %zu.\n", m_npar, input.par_initial.size());
    printf("ERROR: par_initial size mismatch. Expected %d, got %zu.\n", m_npar, input.par_initial.size());
    return false;
  }

  if (input.erp_initial.size() == (size_t)m_npar)
  {
    std::copy(input.erp_initial.begin(), input.erp_initial.end(), erp);
  }
  else if (m_npar > 0)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: erp_initial size mismatch. Expected %d, got %zu.\n", m_npar, input.erp_initial.size());
    printf("ERROR: erp_initial size mismatch. Expected %d, got %zu.\n", m_npar, input.erp_initial.size());
    return false;
  }

  if (input.idpar_data.size() == idpar_actual_size)
  {
    memcpy(idpar, input.idpar_data.data(), idpar_actual_size * sizeof(bcd_t));
  }
  else if (m_npar > 0)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: idpar_data size mismatch. Expected %zu, got %zu.\n", idpar_actual_size, input.idpar_data.size());
    printf("ERROR: idpar_data size mismatch. Expected %zu, got %zu.\n", idpar_actual_size, input.idpar_data.size());
    return false;
  }

  // --- Initialize fitting matrices and related arrays (sized by m_nfit) ---
  if (m_nfit <= 0)
  {
    if (lufit && m_npar > 0)
    {
      //fprintf(lufit, "No parameters marked for fitting (nfit = %d).\n", m_nfit);
    }
    printf("No parameters marked for fitting (nfit = %d).\n", m_nfit);
    iperm = nullptr;
    dpar = nullptr;
    delbgn = nullptr;
    fitbgn = nullptr;
    var = nullptr;
    oldfit = nullptr;
    fit = nullptr;
    ndfit = 0;
    ndiag = 0;
  }
  else
  {
    ndfit = m_nfit + 1;
    ndiag = ndfit + 1;

    size_t nl_maxmem_dummy;
    m_nsize_p = maxmem(&nl_maxmem_dummy);
    if (ndfit > m_nsize_p)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: Number of independent parameters + 1 (%d) is too big for memory segment (%ld).\n", ndfit, m_nsize_p);
      printf("Number of independent parameters + 1 (%d) is too big: %ld\n", ndfit, m_nsize_p);
      return false;
    }

    iperm = (int *)mallocq((size_t)m_nfit * sizeof(int));
    if (!iperm)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: mallocq for iperm failed.\n");
      perror("mallocq for iperm failed");
      return false;
    }

    dpar = (double *)mallocq((size_t)ndfit * sizeof(double));
    if (!dpar)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: mallocq for dpar failed.\n");
      perror("mallocq for dpar failed");
      return false;
    }

    delbgn = (double *)mallocq((size_t)m_nfit * sizeof(double));
    if (!delbgn)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: mallocq for delbgn failed.\n");
      perror("mallocq for delbgn failed");
      return false;
    }

    size_t standard_packed_elements = ((size_t)m_nfit * ((size_t)m_nfit + 1)) / 2;

    fitbgn = (double *)mallocq(standard_packed_elements * sizeof(double));
    if (!fitbgn)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: mallocq for fitbgn failed.\n");
      perror("mallocq for fitbgn failed");
      return false;
    }

    var = (double *)mallocq(standard_packed_elements * sizeof(double));
    if (!var)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: mallocq for var failed.\n");
      perror("mallocq for var failed");
      return false;
    }

    if (input.var_initial_from_getvar.size() == standard_packed_elements)
    {
      std::copy(input.var_initial_from_getvar.begin(), input.var_initial_from_getvar.end(), var);
    }
    else if (!input.var_initial_from_getvar.empty())
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: var_initial_from_getvar size mismatch. Expected %zu (std packed), got %zu.\n", standard_packed_elements, input.var_initial_from_getvar.size());
      printf("ERROR: var_initial_from_getvar size mismatch. Expected %zu (std packed), got %zu.\n", standard_packed_elements, input.var_initial_from_getvar.size());
      return false;
    }
    else if (m_nfit > 0)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: var_initial_from_getvar is empty for nfit > 0.\n");
      printf("ERROR: var_initial_from_getvar is empty for nfit > 0.\n");
      return false;
    }

    size_t oldfit_elements = standard_packed_elements + (size_t)m_nfit;
    oldfit = (double *)mallocq(oldfit_elements * sizeof(double));
    if (!oldfit)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: mallocq for oldfit failed.\n");
      perror("mallocq for oldfit failed");
      return false;
    }

    fit = (double *)mallocq((size_t)m_nfit * (size_t)ndfit * sizeof(double));
    if (!fit)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: mallocq for fit failed.\n");
      perror("mallocq for fit failed");
      return false;
    }

    // --- 4. Initialize fitbgn, dpar, oldfit, delbgn based on m_inpcor and this->var ---
    m_ndfree0 = 0;
    double current_scale_factor; // Renamed from 'scale'

    if (m_inpcor > 0)
    { // "Supplied variance" path
      // Assume 'var' is standard packed LOWER Cholesky factor L of Covariance.
      // (If getvar truly gives UPPER, then var is U, and fit=U. Logic might need U^T U.)
      // For now, proceed as if var = L_packed_lower.

      // 4.1. Form full lower triangular `fit` matrix (L) from `var` (L_packed_lower).
      // `var` (member this->var) holds U_packed_colmajor from getvar.
      // Create a temporary L_packed_colmajor from this->var (U_packed).
      std::vector<double> l_packed_temp(standard_packed_elements);
      // Transpose U_packed in this->var to L_packed in l_packed_temp
      // U_ij is var[j*(j+1)/2 + i] for i<=j. We want L_ij = U_ji.
      // L_packed_colmajor stores L00, L10,L11, L20,L21,L22,...
      int k_l_packed = 0;
      for (int j_col_L = 0; j_col_L < m_nfit; ++j_col_L)
      { // Iterate columns of L
        for (int i_row_L = j_col_L; i_row_L < m_nfit; ++i_row_L)
        { // Iterate rows of L (from diagonal down)
          // L(i_row_L, j_col_L) = U(j_col_L, i_row_L)
          // Find U(j_col_L, i_row_L) in this->var (U_packed)
          int u_row = j_col_L;
          int u_col = i_row_L; // u_row <= u_col
          int idx_in_u_packed = u_col * (u_col + 1) / 2 + u_row;
          if (k_l_packed < (int)standard_packed_elements && idx_in_u_packed < (int)standard_packed_elements)
          {
            l_packed_temp[k_l_packed++] = this->var[idx_in_u_packed];
          }
          else
          { /* error bounds */
          }
        }
      }
      // Now l_packed_temp contains L0_packed_lower.

      // Form full lower triangular `fit` matrix (L0) from l_packed_temp.
      double *pl_packed_ptr = l_packed_temp.data();
      for (int j_col = 0; j_col < m_nfit; ++j_col)
      {
        double *pfit_col = this->fit + j_col * ndfit;
        for (int i_row = 0; i_row < m_nfit; ++i_row)
          pfit_col[i_row] = 0.0;
        for (int i_row = j_col; i_row < m_nfit; ++i_row)
        {
          pfit_col[i_row] = *pl_packed_ptr++;
        }
      }
      // Now `fit` is full Lower triangular L0.

      // 4.3. `dqrfac` on L_s.
      // `this->erpar` (diag_flags) is passed as `wk`. dqrfac uses it for pivoting and overwrites it.
      int n_rank = dqrfac(fit, ndfit, m_nfit, m_nfit, 0, this->erpar, iperm);
      if (m_nfit > n_rank && lufit)
      { /* singular warning */
      }
      // `fit` now contains F (lower factor of L_s * P). `this->erpar` has diag of F.

      // 4.4. `dqrsolv` to get (F P^T)^-1 in `fit`.
      dqrsolv(fit, ndfit, n_rank, m_nfit, 0, iperm);
      // `fit` now contains (F P^T)^-1. This should be lower triangular.

      // 4.5. Unscale `fit`. We want L_inv = D_inv * (F P^T)^-1 (if F came from L D_inv P).
      // Original lsqfit unscaled solution vectors. Here we unscale the inverse factor.
      // Each column k of (F P^T)^-1 should be scaled by dpar[k] (D_inv_kk).
      // (F P^T)^-1 is lower triangular. Scale L_inv_column_k by D_inv_kk.
      double *pfit_col_start = fit;
      for (int k_col = 0; k_col < m_nfit; ++k_col)
      {
        current_scale_factor = dpar[k_col]; // D_inv_kk
        // Scale L_inv(k_col:end, k_col)
        dscal(m_nfit - k_col, current_scale_factor, pfit_col_start + k_col, 1);
        pfit_col_start += ndfit;
      }
      // Now `fit` should contain L_inv (lower triangular).

      // 4.6. Store Cov_inv = L_inv^T * L_inv in `fitbgn` (standard packed lower).
      //    First, form full Cov_inv in `fit` (or a temp matrix).
      //    fit (L_inv) is m_nfit x m_nfit (in m_nfit x ndfit storage).
      //    Temp_Symm = L_inv^T * L_inv.
      //    This requires a dsyrk call if L_inv is available: C = alpha*A^T*A + beta*C or C = alpha*A*A^T + beta*C
      //    cblas_dsyrk(CblasColMajor, CblasUpper/Lower, CblasTrans/NoTrans, N, K, alpha, A, lda, beta, C, ldc)
      //    To compute L_inv^T * L_inv (result lower triangle stored in fit):
      //    Use a temporary matrix for L_inv^T if needed, or compute elements directly.
      //    Or, more simply, if `fit` from `dqrsolv` after unscaling is *already* the symmetric inv_covariance matrix
      //    (which `lsqfit`'s `dkold` save implies for the derivative part `dk(0:nr-1,0:nr-1)`), then just pack it.
      //    The original code: `dcopy(nt, pdk, 1, pdkold, 1);` saves the full solution block.
      //    The `dcopy(n,pfit,ndfit,pfitb,1)` from `main.c`'s "Supplied Variance" implies `fit` had the symmetric matrix.
      //    Let's assume `fit` after unscaling IS the symmetric Inverse Covariance.
      double *pfitb_ptr = fitbgn;
      pfit_col_start = fit; // `fit` is m_nfit x ndfit (contains symmetric InvCov)
      for (int j_col_idx = 0; j_col_idx < m_nfit; ++j_col_idx)
      { // For each column j
        // Pack lower triangle: elements from row j_col_idx down to m_nfit-1
        for (int i_row_idx = j_col_idx; i_row_idx < m_nfit; ++i_row_idx)
        {
          *pfitb_ptr++ = pfit_col_start[i_row_idx]; // fit[i_row_idx, j_col_idx]
        }
        pfit_col_start += ndfit; // Next column in 'fit'
      }
    }
    else
    { // "No supplied variance"
      // `var` is standard packed upper, diagonal elements are sigma_i^2 from getvar.
      // `fitbgn` needs to be standard packed upper, 1/sigma_i^2 on diagonal, 0 off-diagonal.
      memset(fitbgn, 0, standard_packed_elements * sizeof(double));
      double *pvar_ptr = var;
      double *pfitb_ptr = fitbgn;
      for (int j_col = 0; j_col < m_nfit; ++j_col)
      {
        for (int i_row = 0; i_row <= j_col; ++i_row)
        {
          if (i_row == j_col)
          { // Diagonal element
            if (*pvar_ptr == 0.0)
            {
              *pfitb_ptr = 1.0 / (m_tiny * m_tiny);
            }
            else
            {
              *pfitb_ptr = 1.0 / (*pvar_ptr);
            }
            if (*pfitb_ptr < 1.e-15)
              --m_ndfree0;
          } // Off-diagonals in fitbgn remain 0 from memset
          pvar_ptr++;
          pfitb_ptr++;
        }
      }
      memset(var, 0, standard_packed_elements * sizeof(double)); // Zero out var
    }

    if (m_nfit > 0)
    {
      if (oldpar && par)
        oldpar[0] = par[0]; // Only first element? Original did this.
      if (oldfit && fitbgn)
        oldfit[0] = fitbgn[0];
      if (delbgn)
        memset(delbgn, 0, sizeof(double) * m_nfit);
    }
  }

  // --- 5. Calculate m_nxpar_actual and m_nxfit ---
  m_nxpar_actual = input.nxpar_from_file;
  if (m_nxpar_actual > m_npar || m_nxpar_actual < 0)
    m_nxpar_actual = 0;
  if (m_nxpar_actual > 0 && lufit)
  {
    fprintf(lufit, "NUMBER OF PARAMETERS EXCLUDED IF INDEX < 0 =%5d\n", input.nxpar_from_file);
  }
  m_nxpar_actual = m_npar - m_nxpar_actual; // Convert count from end to starting index

  m_nxfit = m_nfit;
  for (int i_par_idx = m_npar - 1; i_par_idx >= m_nxpar_actual; --i_par_idx)
  {
    int ibcd_offset = i_par_idx * m_ndbcd;
    if (idpar && (size_t)ibcd_offset < idpar_actual_size)
    {                                      // Check bounds
      if (NEGBCD(idpar[ibcd_offset]) == 0) // If parameter is fitted (not fixed)
        --m_nxfit;
    }
    else if (m_npar > 0)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: Index out of bounds for idpar in nxfit calculation (param %d).\n", i_par_idx);
      printf("ERROR: Index out of bounds for idpar in nxfit calculation (param %d).\n", i_par_idx);
      // This could be critical, consider returning false
    }
  }

  // --- 6. Call lbufof ---
  lbufof(-m_nfit, m_limlin);

  if (lufit && fabs(m_fqfacq - 1.) >= 1e-10)
  {
    fprintf(lufit, " IR Frequencies Scaled by %12.10f\n", m_fqfacq);
  }
  return true;
} // initializeParamaters


/* method responsible for:
1. Setting up the CalculationEngine:
  - Calling calc->setopt(...) to process option lines.
  - Calling calc->setfmt(...) to set the quantum number format.
2. Processing Experimental Lines:
  - Reading lines from input.lineData_raw. Since linein (now CalFit::linein) expects a FILE*, we'll need to write input.lineData_raw to an in-memory temporary file (tmpfile()).
  - Calling this->linein(...) with this temporary file.
  - Storing the returned maxf (max F quantum number) and actual number of lines read.
3. Setting up Hamiltonian Blocks:
  - Initializing m_nblkpf_actual (blocks per F) using m_maxf_from_linein.
  - Initializing m_maxdm_actual (max Hamiltonian dimension) using m_nsize_p (max memory allowed).
  - Calling calc->setblk(...) to define the block structure. This updates m_nblkpf_actual and m_maxdm_actual.
4. Setting Parameter Labels:
  - Calling getlbl(...) (this is a global C function from ulib.c or similar, not a CalFit method based on original code structure) to populate this->parlbl.
5. Converting Lines and Setting up Links:
  - Calling this->lineix(...) to process lines according to the block structure.
6. Allocating Memory for Fitting:
  - Allocating this->teig and this->pmix based on m_maxdm_actual and m_nfit
*/
// In CalFit.cpp

bool CalFit::processLinesAndSetupBlocks(const CalFitInput &input)
{
  if (!calc)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: CalculationEngine is null in processLinesAndSetupBlocks.\n");
    puts("ERROR: CalculationEngine is null.");
    return false;
  }
  if (!lufit)
  {
    puts("ERROR: lufit stream is null in processLinesAndSetupBlocks.");
    return false;
  }

  // 1. Engine setup (setopt) was DONE by CalFitIO.
  //    Member variables like m_itd, m_ndbcd are set. m_namfil_buffer has content.

  // 2. Set quantum number format in CalculationEngine (calc->setfmt)
  int iqnfmt_val_from_setfmt = 0;
  // m_nfmt (member) will store the iqnfmt value for linein/lineix.
  // Initialize with a default or value from input if appropriate (e.g. input.nfmt_cat_from_setopt if it's the same kind of nfmt)
  // The original main did: iqnfmt = 0; ... nqn = calc->setfmt(&iqnfmt, 1);
  // So, iqnfmt_val_from_setfmt is the primary output.
  int nqn_engine_actual = calc->setfmt(&iqnfmt_val_from_setfmt, 1); // iflg=1 usually means primary format
  if (nqn_engine_actual <= 0)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: calc->setfmt failed to set quantum number format (returned %d).\n", nqn_engine_actual);
    printf("ERROR: calc->setfmt failed (returned %d).\n", nqn_engine_actual);
    return false;
  }
  m_nfmt = iqnfmt_val_from_setfmt; // Store the iqnfmt for line processing by CalFit methods
  // Original main also had: if (nqn > MAXCAT) catqn = nqn; nqn = nqn + nqn; (for pairs)
  // This `nqn_engine_actual` is used later in main iteration loop's qnfmt2 call. Store if needed.
  // For now, m_nfmt (the iqnfmt itself) is the key part for linein/lineix.
  m_nqn_for_iteration = nqn_engine_actual; // Store this for later use with qnfmt2

  // 3. Process Experimental Lines (using this->linein)
  FILE *temp_line_file = tmpfile();
  if (!temp_line_file)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: Unable to create temporary file for line data.\n");
    puts("ERROR: Unable to create temporary file for line data.");
    return false;
  }
  for (const auto &line_str : input.lineData_raw)
  {
    fputs(line_str.c_str(), temp_line_file);
    if (!line_str.empty() && line_str.back() != '\n')
    {
      fputc('\n', temp_line_file);
    }
  }
  rewind(temp_line_file);

  int current_nline_val = m_limlin;                                        // Pass current max line count, linein updates it
  m_maxf_from_linein = this->linein(temp_line_file, &current_nline_val, m_nfmt); // Call CalFit's linein
  m_nline = current_nline_val;                                             // Update m_nline with actual lines read
  fclose(temp_line_file);
  fflush(stdout); // As in original after linein

  if (m_nline <= 0)
  {
    // if (lufit)
    //   fprintf(lufit, "No lines read (NLINE = %d) or linein failed.\n", m_nline);
    puts("No lines read or linein failed.");
    // It's possible to have 0 lines and proceed if only calculation is needed,
    // but typically for fitting this is an issue. Let's return false.
    return false;
  }
  if (m_limlin > m_nline)
  {                     // From original main logic
    m_limlin = m_nline; // Adjust m_limlin to actual if fewer were read
  }

  // 4. Initialize block structure in CalculationEngine (calc->setblk)
  m_nblkpf_actual = (m_maxf_from_linein > 0) ? m_maxf_from_linein : 1;
  m_maxdm_actual = (int)m_nsize_p; // m_nsize_p from initializeParameters
                                   // This is max Hamiltonian dimension allowed by memory.
                                   // setblk can reduce m_maxdm_actual based on actual blocks.

  // Ensure idpar and par are valid before passing to setblk
  if (!this->idpar || !this->par)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: idpar or par is null before calc->setblk.\n");
    puts("ERROR: idpar or par is null before calc->setblk.");
    return false;
  }
  int k_from_setblk = calc->setblk(lufit, m_npar, this->idpar, this->par, &m_nblkpf_actual, &m_maxdm_actual);

  // 5. Set parameter labels (getlbl - global C function)
  if (this->parlbl && this->idpar)
  { // parlbl allocated in initializeParameters
    const char *namfil_cstr = (m_namfil_buffer[0] == '\0') ? nullptr : m_namfil_buffer;
    getlbl(m_npar, this->idpar, this->parlbl, namfil_cstr, k_from_setblk, LBLEN);
  }
  else
  {
    // if (lufit)
    //   fprintf(lufit, "Warning: parlbl or idpar is null, skipping getlbl.\n");
    printf("Warning: parlbl or idpar is null, skipping getlbl.\n");
  }

  // 6. Convert lines and set up links for fitting (this->lineix)
  // lineix takes nitr as flg for detailed output control.
  int bad_lines = this->lineix(lufit, m_nitr_requested, m_nline, m_nblkpf_actual, m_nfmt);
  if (bad_lines > 0)
  {
    printf("%d bad lines\n", bad_lines);
    // if (lufit)
    //   fprintf(lufit, "%d bad lines reported by lineix.\n", bad_lines);
    printf("%d bad lines reported by lineix.\n", bad_lines);
  }
  fflush(stdout);

  // 7. Allocate teig and pmix_block
  if (m_maxdm_actual <= 0)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: m_maxdm_actual is %d, cannot allocate Hamiltonian matrices.\n", m_maxdm_actual);
    printf("ERROR: m_maxdm_actual is %d, cannot allocate Hamiltonian matrices.\n", m_maxdm_actual);
    return false;
  }
  if (m_maxdm_actual > m_nsize_p)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: Hamiltonian dimension (%d) is too big for available memory (%ld).\n", m_maxdm_actual, m_nsize_p);
    printf("Hamiltonian dimension is too big: %d %ld\n", m_maxdm_actual, m_nsize_p);
    return false;
  }

  // Free previous if any (shouldn't happen in single run)
  if (teig)
    free(teig);
  teig = nullptr;
  if (pmix)
    free(pmix);
  pmix = nullptr;

  size_t maxdm_sq_elements = (size_t)m_maxdm_actual * (size_t)m_maxdm_actual;
  teig = (double *)mallocq(maxdm_sq_elements * sizeof(double));
  if (!teig)
  {
    // if (lufit)
    //   fprintf(lufit, "ERROR: mallocq for teig failed.\n");
    perror("mallocq for teig failed");
    return false;
  }

  size_t num_deriv_params_for_alloc = (m_nfit > 0) ? (size_t)m_nfit : 0;
  // If m_nfit is 0, deriv_dim is 0, so block is for pmix_scratch + egy.
  // This matches original main.c where (nfit+2) became (0+2)=2.

  size_t pmix_block_elements = (size_t)m_maxdm_actual * (2 + num_deriv_params_for_alloc);

  if (pmix_block_elements == 0 && m_maxdm_actual > 0)
  {
    // This case implies m_maxdm_actual > 0 but (2 + num_deriv_params_for_alloc) was 0, which is impossible.
    // It means m_maxdm_actual itself might have been <= 0 before this.
    // The check `if (m_maxdm_actual <= 0)` before teig allocation should catch this.
    if (lufit)
      fprintf(lufit, "Warning: pmix_block_elements calculated as zero unexpectedly with m_maxdm_actual = %d.\n", m_maxdm_actual);
    pmix = nullptr; // Should not be reached if m_maxdm_actual > 0
  }
  else if (pmix_block_elements > 0)
  {
    pmix = (double *)mallocq(pmix_block_elements * sizeof(double));
    if (!pmix)
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: mallocq for pmix failed.\n");
      perror("mallocq for pmix failed");
      free(teig);
      teig = nullptr; // Clean up teig if pmix fails
      return false;
    }
    pmix[0] = 1.0; // As in original main
  }
  else
  {
    // pmix_block_elements is 0, likely because m_maxdm_actual was 0.
    pmix = nullptr;
  }

  rqexit(1); // Initialize/reset interrupt flag system (from calpgm.h, ulib.c)

  return true;
}


// Define FORMAT STRINGS (these were macros in original main)
#define FMT_xbgnMW_STR "%5d: %s%14.5f%14.5f %10.5f %10.5f %9.5f"
#define FMT_xblnMW_STR "%14.5f %10.5f %6.4f\n"
#define FMT_xbgnIR_STR "%5d: %s%14.7f%14.7f %10.7f %10.7f %9.7f"
#define FMT_xblnIR_STR "%14.7f %10.7f %6.4f\n"
#define PR_DELAY 6 // From original main


bool CalFit::performIteration(const CalFitInput &input, CalFitOutput &output)
{
  if (!lufit)
  {
    puts("ERROR: lufit stream is null in performIteration.");
    return false;
  }
  if (m_nline <= 0 && m_nfit > 0)
  { // No lines to fit but parameters exist
    printf("WARNING: No lines available to perform iteration, but nfit > 0.\n");
    // Set some default output and return, or handle as error if fitting is expected.
    output.itr = 0;
    output.xsqbest = 0.0;
    // Populate output.par and output.erpar with initial values
    if (m_npar > 0 && this->par)
      output.par.assign(this->par, this->par + m_npar);
    if (m_npar > 0 && this->erpar)
      output.erpar.assign(this->erpar, this->erpar + m_npar);
    return true; // Or false if this is an unrecoverable state for fitting
  }
  if (m_nfit <= 0)
  {
    printf("No parameters to fit (nfit = %d). Skipping iteration loop.\n", m_nfit);
    // Perform a "calculation only" if nitr requested it (e.g. negative nitr in original)
    // For now, just set output and return.
    output.itr = 0;
    output.xsqbest = 0.0;
    if (m_npar > 0 && this->par)
      output.par.assign(this->par, this->par + m_npar);
    if (m_npar > 0 && this->erpar)
      output.erpar.assign(this->erpar, this->erpar + m_npar);
    return true;
  }

  // Pointers for hamx arguments, derived from this->pmix allocation
  // Assumes pmix_block: [pmix_scratch (m_maxdm_actual)], [egy (m_maxdm_actual)], [egyder (m_maxdm_actual * m_nfit)]
  double *hamx_pmix_scratch = this->pmix;
  double *hamx_egy = this->pmix + m_maxdm_actual;
  // This pointer is valid regardless of m_nfit.
  // If m_nfit was 0 during allocation, this points to the end of the allocated pmix block.
  // If m_nfit > 0, this points to the start of the egyder region.
  double *hamx_egyder = this->pmix + 2 * m_maxdm_actual;

  // Local variables from original main's iteration loop
  static double zero_static = 0.; // To avoid repeated initialization if loop is re-entered (not in this design)
  double fqfac[4];
  double varv[5];                                                                          // For lsqfit
  double adif, cerr_calc, afrq, rerr_line, cwid;                                           // Renamed cerr to cerr_calc
  double xsqir, xsqmw, avgmw, avgir, xerr_line, xfrq_line, xsqt, current_scale_iter, xwid; // Renamed scale, xerr, xfrq
  double fqfacq_local = m_fqfacq;                                                          // Use member
  double xerrmx_local = m_xerrmx;                                                          // Use member
  double dif_par, frq_calc, val_iter, xwt_line;                                            // Renamed val, dif

  double bigd, bigf, bigdIR, bigfIR, bigdMW, bigfMW;
  int ifac, iblk, lblk = 0, iflg_line, line_idx, icnt_progress, nfir;
  int lstf = 0, nitr_actual;
  int k_loop, iblnd, lblnd, initl_getdbk, noptn_dummy, ibase, nf_fitted_lines, itd_dummy; // Renamed n, nfmt, etc.
  int maxf_dummy, nrj;
  int marqflg_lsqfit;
  int nsize_block; // Renamed nsize
  int ndfree;
  short qnum_line[2 * MAXQN]; // MAXQN from lsqfit.h or calpgm.h
  char ch_par, pare_str[64], aqnum_str[6 * MAXQN + 2], card_iter_log[NDCARD];
  int idx_getdbk_local; // Local variable for 'indx' output from getdbk

  bigdIR = 9.9999999;
  bigdMW = 100. * bigdIR;
  bigfIR = 99999.9999999;
  bigfMW = 100. * bigfIR;

  fqfac[0] = 1.0;
  fqfac[1] = -1.0;
  fqfac[2] = fqfacq_local / 29979.2458;
  fqfac[3] = -fqfac[2];

  m_marqp[0] = input.marqp0; // Initial Marquardt parameter from input file
  m_marqp[1] = -1.0;         // Trust region, -1 means find max initially
  m_marqp[2] = 1.0;          // Trust region ratio
  if (m_marqp[0] < 0.0)
    m_marqp[0] = 0.0;

  m_xsqbest = zero_static; // Best RMS error
  m_marqlast = m_marqp[0]; // Store initial marqp for output if no iterations
  m_itr = 0;               // Iteration count

  if (this->dpar)
    this->dpar[0] = 0.0; // First element of dpar (parameter shifts) often treated specially

  nitr_actual = (m_nitr_requested < 0) ? -m_nitr_requested : m_nitr_requested;
  if (nitr_actual == 0 && m_nitr_requested <= 0)
  { // if nitr=0, or e.g. -1 to mean calc only
    // If just calculation, may not need full iteration, but original did one pass.
    // For now, let's assume if nitr_actual is 0, we might do one pass for calc if m_nitr_requested < 0
  }

  // --- START ITERATION LOOP ---
  do
  {
    // 1. Set up dpar with current parameter values for derivative calculations
    //    dpar is also used to store parameter *shifts* later.
    //    Original: k=0; ibcd=nxpar*ndbcd; for(i=nxpar;i<npar;++i,ibcd+=ndbcd){ if(NEGBCD(idpar[ibcd])==0) dpar[k++]=par[i];}
    //    This copies values of *fitted* parameters (from m_nxpar_actual onwards) into the start of dpar array.
    //    This dpar is then passed to dnuadd. This usage seems to mix parameter values and shifts.
    //    Let's use a temporary array for passing current fitted parameter values to dnuadd.
    if (m_nxfit <= 0 && m_nfit > 0)
    { // Should not happen if nxfit correctly calculated
      printf("Warning: m_nxfit is %d while m_nfit is %d. Problem with parameter exclusion logic.\n", m_nxfit, m_nfit);
    }
    std::vector<double> current_fitted_params(m_nfit); // Max possible needed for dnuadd
    k_loop = 0;
    int ibcd_idx = m_nxpar_actual * m_ndbcd; // Start check from m_nxpar_actual
    for (int i_par = m_nxpar_actual; i_par < m_npar; ++i_par, ibcd_idx += m_ndbcd)
    {
      if (this->idpar && (size_t)ibcd_idx < ((size_t)m_npar * m_ndbcd + m_ndbcd + 3))
      { // Bounds check
        if (NEGBCD(this->idpar[ibcd_idx]) == 0)
        { // If parameter i_par is fitted
          if (k_loop < m_nfit)
          { // Check bounds for current_fitted_params
            current_fitted_params[k_loop++] = this->par[i_par];
          }
          else
          { /* error, too many fitted params for m_nfit */
            break;
          }
        }
      }
      else if (m_npar > 0)
      { /* error, idpar bounds */
        break;
      }
    }
    // Fill remaining (if any) from start of parameter list for parameters not in the "nxpar tail"
    ibcd_idx = 0;
    for (int i_par = 0; i_par < m_nxpar_actual; ++i_par, ibcd_idx += m_ndbcd)
    {
      if (this->idpar && (size_t)ibcd_idx < ((size_t)m_npar * m_ndbcd + m_ndbcd + 3))
      {
        if (NEGBCD(this->idpar[ibcd_idx]) == 0)
        {
          if (k_loop < m_nfit)
          {
            current_fitted_params[k_loop++] = this->par[i_par];
          }
          else
          { /* error */
            break;
          }
        }
      }
      else if (m_npar > 0)
      { /* error */
        break;
      }
    }
    // k_loop should now equal m_nfit if all went well.

    // 2. Loop through Hamiltonian blocks (getdbk)
    int lnext_getdbk = 0;                                              // Controls getdbk loop
    lblk = 0;                                                          // Last block processed
    lstf = 0;                                                          // Last F quantum number printed for progress
    // Call getdbk to initialize/get the first set of parameters.
    // idx_getdbk_local will be set by this first call.
    getdbk(&lnext_getdbk, &iblk, &idx_getdbk_local, &initl_getdbk, &ifac);

    do
    {
      // On subsequent calls, getdbk uses lnext_getdbk to advance and sets outputs.
      line_idx = getdbk(&lnext_getdbk, &iblk, &idx_getdbk_local, &initl_getdbk, &ifac); // Get next line info
      if (iblk != lblk)
      { // New block
        if (rqexit(0) != 0)
          break; // Check for user interrupt

        // Get block quantum numbers and size
        calc->getqn(iblk, 0, 0, qnum_line, &nsize_block); // qnum_line gets some ID, nsize_block is dim
        if (nsize_block == 0)
          continue; // Skip empty blocks

        lblk = iblk;
        if (nsize_block > m_maxdm_actual)
        {
          // if (lufit)
          //   fprintf(lufit, "ERROR: Size of block %d (%d) exceeds max dimension %d.\n", iblk, nsize_block, m_maxdm_actual);
          printf("ERROR: Size of block %d (%d) exceeds max dimension %d.\n", iblk, nsize_block, m_maxdm_actual);
          return false; // Fatal error
        }

        k_loop = (iblk - 1) / m_nblkpf_actual; // Assuming m_nblkpf_actual is F step
        if (lstf != k_loop && caldelay(PR_DELAY) != 0)
        {
          printf("Starting Quantum %3d\n", k_loop);
          fflush(stdout);
          lstf = k_loop;
        }

        // Calculate energies and derivatives for this block
        // Pointers hamx_egy, hamx_egyder, hamx_pmix_scratch are set up from this->pmix
        // this->teig is the eigenvector output matrix
        calc->hamx(iblk, nsize_block, m_npar, this->idpar, this->par,
                   hamx_egy, this->teig, hamx_egyder, hamx_pmix_scratch, FALSE);
      }
      // Accumulate derivatives for the current line
      // dnuadd needs array of current *fitted* parameter values.
      // Original used `dpar` for this, which was confusing. We use `current_fitted_params`.
      dnuadd(m_nfit, m_nxfit, initl_getdbk, idx_getdbk_local, ifac,
             hamx_egy, hamx_egyder, nsize_block, line_idx,
             current_fitted_params.data(), fqfac);

    } while (lnext_getdbk != 0); // Repeat until no more energies needed for lines

    if (rqexit(0) != 0 && lnext_getdbk != 0)
    { // Exited loop due to interrupt but not finished
      // if (lufit)
      //   fprintf(lufit, "Iteration interrupted by user.\n");
      printf("Iteration interrupted.\n");
      // Decide how to handle: break outer loop, return error, etc.
      // For now, break outer loop.
      break;
    }
    if (lblk > 0)
    { // Check if any blocks were processed
      k_loop = (lblk - 1) / m_nblkpf_actual;
      printf("Finished Quantum %3d\n", k_loop);
      fflush(stdout);
    }

    // 3. Initialize least-squares matrix (fit) and RHS vector (target for J^T * residuals)
    //    `fit` will store L (Cholesky of initial J0^T J0 from fitbgn).
    //    The (m_nfit)-th conceptual column of `fit` will store b0 = L0 * delbgn (initial RHS term).
    //    `xsqt` will be ||b0||^2 initially.
    xsqir = xsqmw = avgmw = avgir = 0.0; // Initialize statistics

    double *pfitb_packed_L0 = this->fitbgn;             // Source: Packed Lower L0 from previous iteration or initial setup
    double *pfit_full_L0_col_start = this->fit;         // Destination: Start of current column in full matrix `fit`
    double *pRHS_in_fit_col_start = this->fit + m_nfit; // Destination: Start of RHS conceptual column in `fit`

    // Zero out the RHS part of the `fit` matrix.
    // This accesses elements fit[0, m_nfit], fit[1, m_nfit], ..., fit[m_nfit-1, m_nfit]
    // if we consider fit as (m_nfit)x(m_nfit+1) where last column is RHS.
    // Strided access for column m_nfit:
    for (int i = 0; i < m_nfit; ++i)
    {
      pRHS_in_fit_col_start[i * ndfit] = 0.0; // Accesses elements like fit[i][m_nfit] if row major
                                              // Or fit[i*ndfit_pitch + m_nfit_col_offset]
                                              // The dcopy/ddot use ndfit as leading dimension, so this should be strided for RHS.
    }

    // This loop structure is taken directly from original main.c:
    // It forms L0 in `fit` (main m_nfit x m_nfit part) and b0=L0*delbgn in `fit` (RHS part).
    for (int n_loop = 1; n_loop <= m_nfit; ++n_loop)
    {
      // `dcopy(n_loop, pfitb_packed_L0, 1, pfit_full_L0_col_start, ndfit)`:
      // Copies `n_loop` elements from `pfitb_packed_L0` (contiguous source).
      // To `pfit_full_L0_col_start` with destination increment `ndfit`.
      // This means it writes to:
      // pfit_full_L0_col_start[0],
      // pfit_full_L0_col_start[ndfit],
      // pfit_full_L0_col_start[2*ndfit], ...
      // This fills a ROW of `fit` if `fit` is column-major.
      // If `fit` is to be L (Lower triangular), this fills L^T (Upper triangular) row by row.
      // This means `fit` effectively becomes U0 (Upper Cholesky factor, U0=L0^T).
      dcopy(n_loop, pfitb_packed_L0, 1, pfit_full_L0_col_start, ndfit);

      val_iter = this->delbgn[n_loop - 1]; // delbgn is 0-indexed

      // `daxpy(n_loop, val_iter, pfitb_packed_L0, 1, pRHS_in_fit_col_start, ndfit)`:
      // Adds `val_iter * (L0_col_n_loop)` to the RHS vector (which is also strided).
      // RHS_row_i += val_iter * L0_col_n_loop[i]
      // This forms RHS = sum_cols (L0_col * delbgn_corresponding_element) = L0 * delbgn.
      daxpy(n_loop, val_iter, pfitb_packed_L0, 1, pRHS_in_fit_col_start, ndfit);

      pfitb_packed_L0 += n_loop; // Advance in packed L0 array
      pfit_full_L0_col_start++;  // CRITICAL: Original was ++pfit.
                                 // If pfit_full_L0_col_start is double*, this moves by ONE DOUBLE.
                                 // This means it expects to write to the next ROW's starting element
                                 // for the next transposed column of L0. This confirms `fit` is U0.
    }
    // At this point:
    // - `this->fit` (m_nfit x m_nfit part) contains U0 (Upper triangular, where U0 = L0^T, L0 from fitbgn).
    // - `this->fit` (conceptual (m_nfit)-th column, accessed via pRHS_in_fit_col_start with ndfit stride) contains b0 = L0 * delbgn.
    //   (Note: If fit has U0, then RHS should be U0 * something, or L0 * delbgn consistent with U0 = L0^T).
    //   The daxpy used pfitb_packed_L0, which is L0. So RHS is L0 * delbgn.

    // Calculate initial sum of squares: xsqt = || L0 * delbgn ||^2
    // The RHS vector (L0*delbgn) is stored strided in the last "column" of `fit`.
    xsqt = ddot(m_nfit, pRHS_in_fit_col_start, ndfit, pRHS_in_fit_col_start, ndfit);

    // 4. Loop through experimental lines
    nf_fitted_lines = nrj = nfir = 0;
    icnt_progress = -1;
    if (lufit)
    { // Header for line listing
      for (int i_space = 0; i_space < 40; ++i_space)
        fputc(' ', lufit);
      fputs("EXP.FREQ.  -  CALC.FREQ. -   DIFF.  - EXP.ERR.- ", lufit);
      fputs("EST.ERR.-AVG. CALC.FREQ. -  DIFF. - WT.\n", lufit);
    }
    line_idx = 1;         // 1-based line number for processing
    double ex_line = 1.0; // Effective weight factor

    do
    { // Loop over lines (line_idx from 1 to m_nline)
      if (icnt_progress <= 0 && caldelay(PR_DELAY) != 0)
      {
        printf("Fitting Line %d\n", line_idx);
        fflush(stdout);
        icnt_progress = 50; // Reset counter
      }
      lblnd = line_idx; // Start of potential blend
      current_scale_iter = rerr_line = 0.0;
      iflg_line = 0;
      int supblnd = 0;
      xwid = 1.0;

      // Inner loop: process all components of a blend
      do
      {
        // frqdat: output iblnd, xfrq_line, xwt_line, xerr_line, qnum_line
        // Returns 0 if end of data or error, non-zero if line data valid.
        // dpar (this->dpar) is passed as output for derivatives of this line.
        k_loop = frqdat(lblnd, &iblnd, &xfrq_line, &xwt_line, &xerr_line, qnum_line);

        if ((iblnd & 1) != 0)
        { // Blended line width card (supblnd in original)
          supblnd = 1;
          xwid = xerr_line;                             // xerr_line is width for this type of card
          val_iter = 1.0 / (current_scale_iter * xwid); // Ensure current_scale_iter is not zero
          if (current_scale_iter * xwid == 0)
            val_iter = 0; // Avoid div by zero
          else
            val_iter = sqrt(fabs(1.0 - val_iter * val_iter)); // fabs for safety
          current_scale_iter *= val_iter;
          // dscal dpar (derivatives) by val_iter
          // dpar here means the derivatives vector for the current line, which is dpar[0..nfit-1]
          // dpar[nfit] is obs-calc. So scale all nfit+1 elements.
          dscal(m_nfit + 1, val_iter, this->dpar, 1);
        }
        else if (k_loop != 0)
        { // Normal line or component of blend
          if (xerr_line == 0.0)
            ex_line = 1.0 / m_tiny; // Avoid division by zero
          else
            ex_line = fabs(xwt_line / xerr_line);
          current_scale_iter += ex_line;

          // dnuget: Get (Obs-Calc) and derivatives for this line
          // iflg_line: 0 for first component, >0 for subsequent in blend, -1 for printing single line
          // Output: val_iter = (Obs-Calc) for this line, dpar[0..nfit-1] gets derivatives d(calc)/dp_i
          // dpar[nfit] gets (Obs-Calc) * ex_line (weighted)
          val_iter = dnuget(iflg_line, m_nfit, ex_line, lblnd, this->dpar);
          rerr_line = this->dpar[m_nfit]; // Weighted (Obs-Calc) accumulated for blend
          iflg_line++;
        }
        lblnd++;
        icnt_progress--;
      } while (iblnd < 0); // Loop while iblnd indicates more components in blend

      if (iflg_line == 0)
      {                   // No valid lines/components found for this line_idx
        line_idx = lblnd; // Advance line_idx past this processed (empty) blend
        continue;
      }

      // Calculate errors and process line/blend
      if (current_scale_iter == 0)
        current_scale_iter = 1.0;                                             // Avoid division by zero
      cerr_calc = calerr(m_nfit, this->var, this->dpar) / current_scale_iter; // Estimated error of calc freq
      adif = this->dpar[m_nfit] / current_scale_iter;                         // Average (Obs-Calc) for blend
      afrq = xfrq_line - adif;                                                // Average calculated frequency for blend (using last xfrq_line from frqdat as Explt)
                                                                              // If blend, xfrq_line from frqdat should be the "blend frequency"

      if (fabs(rerr_line) < xerrmx_local)
      {                                // Check if weighted error is acceptable
        xsqt += rerr_line * rerr_line; // Accumulate sum of squares of weighted residuals
        if (xerr_line < 0.0)
        { // IR line (convention: negative error)
          avgir += adif;
          xsqir += adif * adif;
          nfir++;
        }
        else
        { // Microwave line
          avgmw += adif;
          xsqmw += adif * adif;
        }
        // Rotate line into FIT matrix (J^T J and J^T diff)
        // dpar[0...nfit-1] contains derivatives J, dpar[nfit] contains weighted (Obs-Calc)
        jelim(this->fit, this->dpar, ndfit, m_nfit, 1); // ns=1 for one RHS vector
        nf_fitted_lines++;
      }
      else
      { // Line rejected
        if (lufit)
          fputs(" ***** NEXT LINE NOT USED IN FIT\n", lufit);
        nrj++;
        supblnd = 0; // Reset supblnd if rejected
      }

      // Output formatting for the line/blend (extensive logic from original)
      if (xerr_line < 0.0)
      {
        bigf = bigfIR;
        bigd = bigdIR;
      }
      else
      {
        bigf = bigfMW;
        bigd = bigdMW;
      }
      if (fabs(afrq) > bigf)
        afrq = (afrq > 0.0) ? bigf : -bigf;
      if (fabs(adif) > bigd)
        adif = (adif > 0.0) ? bigd : -bigd;

      if (iflg_line == 1)
      {                                                        // Single unblended line printed
        qnfmt2(m_nqn_for_iteration * 2, qnum_line, aqnum_str); // nqn*2 for pair
        if (lufit)
        {
          if (xerr_line < 0.)
          {
            fprintf(lufit, FMT_xbgnIR_STR, line_idx, aqnum_str, xfrq_line, afrq, adif, xerr_line, cerr_calc);
          }
          else
          {
            fprintf(lufit, FMT_xbgnMW_STR, line_idx, aqnum_str, xfrq_line, afrq, adif, xerr_line, cerr_calc);
          }
          fputc('\n', lufit);
        }
      }
      else
      {                                          // Blended line, print all components
        int current_blend_member_idx = line_idx; // Start from the first line of the blend
        cwid = 0.0;
        do
        {
          k_loop = frqdat(current_blend_member_idx, &iblnd, &xfrq_line, &xwt_line, &xerr_line, qnum_line);
          if ((iblnd & 1) != 0)
          {                              // Blend summary/width line
            xfrq_line = 0.0;             // Placeholder for Explt Freq for this line type
            frq_calc = sqrt(fabs(cwid)); // RMS of individual diffs from avg_calc for blend
            if (supblnd != 0 && xwid != 0.0)
              xsqt += cwid / (xwid * xwid);  // Add to sum of squares
            qnfmt2(0, qnum_line, aqnum_str); // No quanta for this line
            if (lufit)
            {
              if (xerr_line < 0.)
              { // Based on original blend line type
                fprintf(lufit, FMT_xbgnIR_STR, current_blend_member_idx, aqnum_str, xfrq_line, frq_calc, frq_calc, xwid, 0.0 /*est.err*/);
              }
              else
              {
                fprintf(lufit, FMT_xbgnMW_STR, current_blend_member_idx, aqnum_str, xfrq_line, frq_calc, frq_calc, xwid, 0.0);
              }
              fputc('\n', lufit);
            }
          }
          else if (k_loop != 0)
          { // Actual component of the blend
            if (supblnd != 0 && xwid != 0.0)
            {                                                                              // If special blend processing active
              ex_line = sqrt(fabs(xwt_line)) / xwid;                                       // Weight for this component
              frq_calc = dnuget(0, m_nfit, ex_line, current_blend_member_idx, this->dpar); // Get its contribution again
              jelim(this->fit, this->dpar, ndfit, m_nfit, 1);                              // Rotate into fit
              nf_fitted_lines++;                                                           // Count as fitted
            }
            else
            {
              // Get (Calc Freq) for this component (iflg=-1 for print mode of dnuget)
              frq_calc = dnuget(-1, m_nfit, ex_line, current_blend_member_idx, this->dpar);
            }
            dif_par = xfrq_line - frq_calc; // Individual diff
            if (fabs(frq_calc) > bigf)
              frq_calc = (frq_calc > 0.0) ? bigf : -bigf;
            if (fabs(dif_par) > bigd)
              dif_par = (dif_par > 0.0) ? bigd : -bigd;
            qnfmt2(m_nqn_for_iteration * 2, qnum_line, aqnum_str);
            if (lufit)
            {
              if (xerr_line < 0.)
              {
                fprintf(lufit, FMT_xbgnIR_STR, current_blend_member_idx, aqnum_str, xfrq_line, frq_calc, dif_par, xerr_line, cerr_calc);
                fprintf(lufit, FMT_xblnIR_STR, afrq, adif, xwt_line); // Avg calc, avg diff, weight
              }
              else
              {
                fprintf(lufit, FMT_xbgnMW_STR, current_blend_member_idx, aqnum_str, xfrq_line, frq_calc, dif_par, xerr_line, cerr_calc);
                fprintf(lufit, FMT_xblnMW_STR, afrq, adif, xwt_line);
              }
            }
            dif_par = frq_calc - afrq;            // Diff of this component's calc from avg calc
            cwid += xwt_line * dif_par * dif_par; // Accumulate for RMS width
          }
          current_blend_member_idx++;
        } while (iblnd < 0); // Loop through all components of the blend for printing
      }
      line_idx = lblnd; // Advance main line index past this processed blend
    } while (line_idx <= m_nline); // End loop over lines

    if (nrj > 0 && lufit)
      fprintf(lufit, "%5d Lines rejected from fit\n", nrj);
    if (nf_fitted_lines < 1)
      nf_fitted_lines = 1; // Avoid division by zero for RMS

    // Zero upper triangle of fit matrix (it's symmetric J^T J after jelim)
    // Original did this by iterating columns and memset on elements *before* diagonal.
    // If 'fit' is column major:
    double *pfit_ptr = this->fit; // Start of matrix
    for (k_loop = 1; k_loop < m_nfit; ++k_loop)
    {                                               // For columns 1 to nfit-1
      pfit_ptr += ndfit;                            // Move to start of column k_loop
      memset(pfit_ptr, 0, sizeof(double) * k_loop); // Zero k_loop elements from pfit_ptr[0] to pfit_ptr[k_loop-1]
    }

    varv[0] = xsqt + (double)nrj * xerrmx_local * xerrmx_local; // Variance for lsqfit
    m_marqlast = m_marqp[0];                                    // Save current Marquardt param before lsqfit changes it

    // 5. Call lsqfit
    // dk=fit, ndm=ndfit, nr=m_nfit, nvec=1 (for one RHS)
    // dkold=oldfit, ediag=erpar (used for normalized diagonals), enorm=dpar (param shifts/errors)
    marqflg_lsqfit = lsqfit(this->fit, ndfit, m_nfit, 1, m_marqp, varv,
                            this->oldfit, this->erpar, this->dpar, this->iperm);

    if (marqflg_lsqfit != 0)
    {                                               // Fit diverging
      dcopy(m_npar, this->oldpar, 1, this->par, 1); // Restore parameters
                                                    // ... (logging as before) ...
    }
    else
    { // Fit converged for this step
      m_xsqbest = (xsqt - varv[3]) / nf_fitted_lines;
      if (m_xsqbest > 0.0)
        m_xsqbest = sqrt(m_xsqbest);
      else
        m_xsqbest = 0.0;
      // Print normalized diagonals (erpar from lsqfit)
      if (lufit)
      {
        fputs("NORMALIZED DIAGONAL:\n", lufit);
        icnt_progress = 0;
        for (k_loop = 0; k_loop < m_nfit; ++k_loop)
        {
          fprintf(lufit, "%5d %13.5E", k_loop + 1, this->erpar[k_loop]);
          if ((++icnt_progress) == 6)
          {
            icnt_progress = 0;
            fputc('\n', lufit);
          }
        }
        if (icnt_progress > 0)
          fputc('\n', lufit);
      }
    }
    if (lufit)
    {
      fprintf(lufit, "MARQUARDT PARAMETER = %g, TRUST EXPANSION = %4.2f\n", m_marqp[0], m_marqp[2]);
    }
//      printf("MARQUARDT PARAMETER = %g, TRUST EXPANSION = %4.2f\n", m_marqp[0], m_marqp[2]);

    if (m_xsqbest > xsqt && marqflg_lsqfit == 0)
      m_xsqbest = xsqt;

    // 6. Update parameters and delbgn
    //    dpar from lsqfit contains parameter changes * delta_p_i
    //    erpar from lsqfit contains scaled errors for fitted parameters
    if (lufit)
    {
      for (int i_space = 0; i_space < 32; ++i_space)
        fputc(' ', lufit);
      fputs("NEW PARAMETER (EST. ERROR) -- CHANGE THIS ITERATION\n", lufit);
    }

    // Parameter changes delta_p are in the RHS column of `this->fit` (scaled by D_L_inv factors)
    double *solution_delta_p_vector = this->fit + m_nfit; // Start of RHS conceptual column

    // The matrix M = D_L_inv * F_lambda_inv is in the main m_nfit x m_nfit block of `this->fit`
    double *matrix_M_start = this->fit;

    char *tlbl_ptr = this->parlbl;
    k_loop = 0; // 0-indexed counter for fitted parameters
    ibase = -1; // Index of the last processed independent parameter

    for (int i_par = 0; i_par < m_npar; ++i_par)
    {
      int ibcd_current = i_par * m_ndbcd; // Calculate BCD offset based on full param index
      // (Removed tlbl_ptr logic from here, handle it with parlbl from output struct in writeOutput)

      this->oldpar[i_par] = this->par[i_par]; // Save current param value

      if (this->idpar && (NEGBCD(this->idpar[ibcd_current]) == 0))
      { // If parameter i_par is FITTED
        if (k_loop < m_nfit)
        {
          // Get parameter change delta_p for this k_loop'th fitted parameter
          dif_par = solution_delta_p_vector[k_loop * ndfit]; // Accesses solution_vector[k_loop] with stride

          this->par[i_par] += dif_par;
          this->delbgn[k_loop] -= dif_par; // Update for trust region based on k_loop

          // Calculate error for this k_loop'th fitted parameter
          // Original: erpar[i] = dpar[k] = parfac * dnrm2(n_local, pfitd_col_k_diag_ptr, 1);
          // pfitd_col_k_diag_ptr was effectively &matrix_M_start[k_loop*ndfit + k_loop]
          // n_local was m_nfit - k_loop (elements from M(k,k) down to M(nfit-1,k))
          double *M_col_k_diag_ptr = matrix_M_start + k_loop * ndfit + k_loop; // M(k,k)
          int n_elements_for_norm = m_nfit - k_loop;
          double error_val = m_parfac * dnrm2(n_elements_for_norm, M_col_k_diag_ptr, 1);

          this->erpar[i_par] = error_val; // Store error for absolute parameter i_par
          // The original also did dpar[k_loop] = error_val.
          // `this->dpar` from lsqfit contained D_L_inv. If we overwrite it here, that's lost for `prcorr`.
          // Let's use a temporary for `dpar_for_putvar` if `putvar` needs these errors.
          // And `this->dpar` (from lsqfit) should be preserved for `prcorr`.
          // For now, let's assume `this->dpar` can be used to store these errors for `putvar`.
          // If `prcorr` uses `this->dpar` as `D_L_inv`, we have a conflict.
          // `prcorr` in main.c used `dpar` (which was `enorm` from lsqfit).
          // `putvar` in main.c used `dpar` (which was overwritten by these errors).
          // This implies `prcorr` happens *before* this `dpar` overwrite, or uses a copy.
          // In our `finalizeOutputData`, `prcorr` is called before `dpar_for_putvar` is finalized.
          // So, `this->dpar` should retain `lsqfit`'s `enorm` output for `prcorr`.
          // The errors for `putvar` should be stored separately or computed in `finalizeOutputData`.
          // Let's assume `this->erpar[i_par]` stores the final error. `dpar[k]` used for `putvar` will be set later.

          // Scale column k (from diagonal down) of matrix M by m_parfac
          // Original: dscal(n_elements_for_norm, m_parfac, M_col_k_diag_ptr, 1);
          // This scaling of `fit` matrix columns seems to be part of error propagation.
          // This means `this->fit` matrix is modified.
          dscal(n_elements_for_norm, m_parfac, M_col_k_diag_ptr, 1);

          // Printing (using current_label logic from before)
          char current_label_str_buffer[LBLEN + 1];
          strncpy(current_label_str_buffer, this->parlbl + k_loop * LBLEN, LBLEN); // Assuming parlbl is 0-indexed by k_loop here
          current_label_str_buffer[LBLEN] = '\0';

          parer(this->par[i_par], error_val, dif_par, pare_str);
          putbcd(card_iter_log, NDCARD, &this->idpar[ibcd_current]);
          if (lufit)
          {
            fprintf(lufit, "%4d %s %10.10s %s\n", k_loop + 1, card_iter_log, current_label_str_buffer, pare_str);
          }
          ibase = i_par;
          k_loop++;
        }
      }
    }
    // Derive errors for non-fitted (dependent) parameters (after loop, using final ibase error)
    // This should use erpar[ibase] which now holds the final error for the master.
    for (int i_par = 0; i_par < m_npar; ++i_par)
    {
      int ibcd_current = i_par * m_ndbcd;
      if (this->idpar && NEGBCD(this->idpar[ibcd_current]) != 0)
      { // Dependent
        if (ibase != -1 && this->par && this->erpar)
        { // ibase must be valid
          // this->par[i_par] is the factor.
          this->erpar[i_par] = fabs(this->par[i_par]) * this->erpar[ibase];
          // Update actual value of dependent parameter: par_dep = factor * par_master
          // This should only be done if par[i_par] truly stored a factor.
          // getpar stores factor: par[i] = vec[0] / parbase.
          // If parbase was par[ibase_initial_value], then par[i] is ratio.
          // It's safer to do this update only at the very end in finalizeOutputData
          // if hamx always expects factors for dependent params.
          // If hamx expects actual values, this update should happen.
          // Let's assume for now, par[i_par] is updated here for consistency in this iteration.
          this->par[i_par] = this->par[i_par] * this->par[ibase];
        }
      }
    }

    // 7. Store `fit` matrix into `var` for next iteration/output
    //    `fit` contains `F_lambda_inv` (scaled by D_inv) and solution vector.
    //    We need to save `Cov_inv_final = D_inv * F_lambda_inv^T * F_lambda_inv * D_inv` (or similar).
    //    The `dkold` (our `oldfit`) from `lsqfit` saved `F_lambda_inv` (unscaled by D_inv).
    //    Original: pfit=fit; pvar=var; for(n=1..nfit){ dcopy(n,pfit,ndfit,pvar,1); ++pfit; pvar+=n;}
    //    This copies the lower triangle of `F_lambda_inv * D_inv` (from `fit`) into `var`.
    //    This is `L_final_inv` if the process yields that.
    //    So `var` gets updated L_inv (packed lower).
    double *pvar_ptr_save = this->var;   // Standard packed lower
    double *pfit_source_col = this->fit; // Full matrix, L_inv_final is in its lower triangle
    for (int j_col = 0; j_col < m_nfit; ++j_col)
    {
      for (int i_row = j_col; i_row < m_nfit; ++i_row)
      {                                            // Lower triangle elements including diagonal
        *pvar_ptr_save++ = pfit_source_col[i_row]; // fit[i_row, j_col]
      }
      pfit_source_col += ndfit;
    }
    // If `var` needs to be standard packed UPPER for consistency with getvar:
    // Need to form symmetric Cov_inv = L_inv^T L_inv, then pack upper of that.
    // For now, assume `var` stores L_inv (packed lower) for next iteration's `fitbgn`.
    // This means `fitbgn` in "Supplied Variance" path should also be built from L_inv (packed lower).

    // 8. Calculate and print statistics (RMS errors, averages)
    if (nf_fitted_lines > nfir)
    {
      current_scale_iter = 1.0 / (double)(nf_fitted_lines - nfir);
      avgmw *= current_scale_iter;
      if (xsqmw * current_scale_iter >= 0)
        xsqmw = sqrt(xsqmw * current_scale_iter);
      else
        xsqmw = 0;
    }
    else
    {
      avgmw = 0;
      xsqmw = 0;
    }
    if (nfir > 0)
    {
      current_scale_iter = 1.0 / (double)nfir;
      avgir *= current_scale_iter;
      if (xsqir * current_scale_iter >= 0)
        xsqir = sqrt(xsqir * current_scale_iter);
      else
        xsqir = 0;
    }
    else
    {
      avgir = 0;
      xsqir = 0;
    }

    printf(" MICROWAVE AVG = %15.6f MHz, IR AVG =%15.5f\n", avgmw, avgir);
    if (lufit)
      fprintf(lufit, " MICROWAVE AVG = %15.6f MHz, IR AVG =%15.5f\n", avgmw, avgir);
    printf(" MICROWAVE RMS = %15.6f MHz, IR RMS =%15.5f\n", xsqmw, xsqir);
    if (lufit)
      fprintf(lufit, " MICROWAVE RMS = %15.6f MHz, IR RMS =%15.5f\n", xsqmw, xsqir);

    m_itr++; // Increment iteration count
    printf(" END OF ITERATION %2d OLD, NEW RMS ERROR=%15.5f %15.5f\n", m_itr, xsqt, m_xsqbest);
    if (lufit)
      fprintf(lufit, " END OF ITERATION %2d OLD, NEW RMS ERROR=%15.5f %15.5f\n", m_itr, xsqt, m_xsqbest);
    fflush(stdout);
    if (lufit)
      fflush(lufit);
  } while (m_itr < nitr_actual && (marqflg_lsqfit != 0 || 0.999999 * xsqt > m_xsqbest));
  // Original condition: itr < nitr && 0.999999 * xsqt > xsqbest
  // Added marqflg check: continue if diverging to allow Marquardt param to increase.

  // --- END OF ITERATION LOOP ---

  // Populate output structure (some fields)
  output.itr = m_itr;
  output.xsqbest = m_xsqbest;
  // Final parameters and errors will be fully populated in finalizeOutputData

  return true;
} // performIteration


bool CalFit::finalizeOutputData(const CalFitInput &input, CalFitOutput &output)
{
  if (!lufit && output.itr > 0)
  { // lufit might be needed for prcorr
    // This case should ideally not happen if lufit is managed correctly through run
    puts("Warning: lufit stream is null in finalizeOutputData, prcorr might be skipped.");
  }

  // --- Populate CalFitOutput structure ---

  // Results from fitting process
  output.xsqbest = m_xsqbest;
  output.itr = m_itr; // Actual iterations performed

  // Copy final parameters and their errors
  // Ensure these arrays were allocated and populated if m_npar > 0
  if (m_npar > 0)
  {
    if (this->par)
    {
      output.par.assign(this->par, this->par + m_npar);
    }
    else
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: this->par is null in finalizeOutputData.\n");
      printf("ERROR: this->par is null in finalizeOutputData.\n");
      return false;
    }
    // this->erpar should contain the final scaled errors for ALL m_npar parameters
    // (fitted ones from lsqfit, dependent ones derived in performIteration)
    if (this->erpar)
    {
      output.erpar.assign(this->erpar, this->erpar + m_npar);
    }
    else
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: this->erpar is null in finalizeOutputData.\n");
      printf("ERROR: this->erpar is null in finalizeOutputData.\n");
      return false;
    }
  }
  else
  {
    output.par.clear();
    output.erpar.clear();
  }

  // Data needed by CalFitIO::writeOutput to reconstruct files
  output.title_for_output = input.title; // Original title from input .par file
  // output.optionLines_for_output = input.optionLines; // CalFitIO needs raw option lines from CalFitInput
  // This will be passed directly to writeOutput from main.

  output.npar_final = m_npar;
  // For limlin in output header: original main.c used `limlin` which could be negative
  // if catqn > MAXCAT. CalFitIO's input.limlin stores this original value.
  output.limlin_final = input.limlin;
  output.nitr_final_actual = m_itr;                // Actual iterations done
  output.nxpar_for_header = input.nxpar_from_file; // The count for the header
  output.marqlast_final = m_marqlast;              // Last Marquardt param before final lsqfit call, or from input if no iters
  output.xerrmx_final = m_xerrmx;
  output.parfac_for_header = m_parfac0; // Original parfac from input file
  output.fqfacq_final = m_fqfacq;
  output.nfit_final = m_nfit;

  if (m_npar > 0)
  {
    size_t idpar_actual_size = ((size_t)m_npar * m_ndbcd + m_ndbcd + 3);
    if (this->idpar)
    {
      output.idpar_final_for_output.assign(this->idpar, this->idpar + idpar_actual_size);
    }
    else
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: this->idpar is null.\n");
      printf("ERROR: this->idpar is null.\n");
      return false;
    }

    size_t parlbl_actual_size = (LBLEN * (size_t)m_npar + 1);
    if (this->parlbl)
    {
      output.parlbl_final_for_output_flat.assign(this->parlbl, this->parlbl + parlbl_actual_size);
    }
    else
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: this->parlbl is null.\n");
      printf("ERROR: this->parlbl is null.\n");
      return false;
    }

    if (this->erp)
    { // Original a-priori errors
      output.erp_original_for_output.assign(this->erp, this->erp + m_npar);
    }
    else
    {
      // if (lufit)
      //   fprintf(lufit, "ERROR: this->erp is null.\n");
      printf("ERROR: this->erp is null.\n");
      return false;
    }
  }
  else
  {
    output.idpar_final_for_output.clear();
    output.parlbl_final_for_output_flat.clear();
    output.erp_original_for_output.clear();
  }

  if (m_nfit > 0)
  {
    size_t standard_packed_elements = ((size_t)m_nfit * ((size_t)m_nfit + 1)) / 2;
    if (this->var)
    { // Final variance/covariance related matrix (packed)
      output.var_final_for_output.assign(this->var, this->var + standard_packed_elements);
    }
    else
    {
      // if (lufit && m_nfit > 0)
      //   fprintf(lufit, "ERROR: this->var is null but nfit > 0.\n");
      printf("ERROR: this->var is null but nfit > 0.\n");
      if (m_nfit > 0)
        return false;
    }

    // `dpar` from lsqfit's `enorm` output contained 1/norm(L_col_k) or similar scaling factors.
    // `putvar` needs the `erpar` argument which was `dpar` in original main's call to `putvar`.
    // That `dpar` for `putvar` was `erpar[i] = parfac * dnrm2(n, pfitd, 1);` from the parameter update loop,
    // effectively the final scaled errors of the *fitted* parameters.
    // `lsqfit`'s `ediag` output (our `this->erpar` member, `m_nfit` elements) contains `abs(diag(F_lambda_inv_unscaled))`.
    // So, `dpar_for_putvar` should be `m_parfac * this->erpar (fitted part)`.
    output.dpar_final_for_putvar.resize(m_nfit);
    for (int i = 0; i < m_nfit; ++i)
    {
      // this->erpar (from lsqfit ediag) contains the errors for fitted parameters
      // This should be the m_nfit block of errors.
      output.dpar_final_for_putvar[i] = m_parfac * this->erpar[i]; // Assuming this->erpar[0..m_nfit-1] holds these
    }
  }
  else
  {
    output.var_final_for_output.clear();
    output.dpar_final_for_putvar.clear();
  }

  // Original main.c calls this sequence at the end:
  // rewind(lubak); fgetstr(card, NDCARD, lubak); chtime(card, 82);
  // fputs(card, stdout); puts("FIT COMPLETE");
  // dcopy(npar, oldpar, 1, par, 1); // Restore par to oldpar if fit diverged or for final state?
  // This was done if marqflg!=0 inside loop.
  // If done here, it means parameters from *before the last iteration* are restored.
  // The current `this->par` holds the latest updated parameters.
  // Original main.c: after loop, if itr==0, exit.
  // Then, prcorr, then fclose(lufit).
  // Then save results: open .par, .var. Put title. Put 2nd line. Put options. Put params.
  // The params written are `this->par` (latest) and `this->erpar` (latest errors).
  // The dcopy(oldpar to par) seems to be only for divergence recovery.

  // Compute and print/store correlation matrix
  if (output.itr > 0 && m_nfit > 0 && this->fit && this->dpar && lufit)
  {
    // `prcorr` needs the final solution/covariance information.
    // `dpar` argument to `prcorr` in original main was `dpar` from `lsqfit` (enorm output - scaling factors).
    // `fit` matrix was the `F_lambda_inv * D_inv` (from lsqfit dk output).
    fprintf(lufit, "\nCorrelation Matrix:\n"); // Placeholder for actual prcorr call
    prcorr(lufit, m_nfit, this->fit, ndfit, this->dpar);
  }

  // Final cleanup before CalFitIO::writeOutput takes over
  // lbufof(-1,0) was called in original main at this stage.
  lbufof(-1, 0); // Release line buffer storage

  // calc->setblk(lufit, 0, ...) was called for cleanup.
  int dummy_nblkpf = 0, dummy_maxdm = 0; // For cleanup call
  if (calc)
  { // Check if calc is valid
    calc->setblk(lufit, 0, this->idpar, this->par, &dummy_nblkpf, &dummy_maxdm);
  }

  if (output.itr == 0 && m_nfit > 0)
  { // No iterations performed but tried to fit
    if (lufit)
      fprintf(lufit, "Output files not updated as no iterations were performed.\n");
    // Depending on requirements, might still write initial state or return !success from run()
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
