/*  supplimentary routines for CALFIT : SUBFIT.C */

/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include "splib/calpgm_types.h"
#include "splib/blas_compat.h"
#include "splib/ulib.h"
#include "spfit/subfit.h"
#include "common/CalError.hpp"
#include "common/file_helpers.hpp"

#define NCARD 130

int getlbl(int npar, bcd_t *idpar, char *parlbl, const char *fil, int idiv, int lblen)
{

  static int ipwr[] = {0, 100, 10000};
  FILE *lu;
  size_t nlbl;
  int kbgn, kend, jdiv, i, k, ibcd, idec, jdec, ibgn, ndbcd, nflg;
  bcd_t *idcmp;
  char buf[NCARD], *parx, *pstr;
  const char *pname;

  /*  gets parameter label from user file */

  /*     IDPAR = list of parameter ID numbers */
  /*     PARLBL = returned list of parameter labels */
  /*     LU = unit to use for user file */
  /*     FIL = file name to use for userfile */
  /*     IDIV = divisor for vibrational field */

  kbgn = -1; kend = 0; nlbl = (size_t) lblen;
  jdiv = 0; jdec = 0; pstr = parlbl;
  for (i = 0; i < npar; ++i) {
    if (*pstr == '\0') {
      kend = i;
      if (kbgn < 0)
        kbgn = i;
    }
    pstr += nlbl;
  }
  if (kbgn < 0)
    return 0;
  lu = NULL;
  pname = fil;  buf[0] = '\0';
  for (idec = 0; idec < 3; ++idec) {
    if (idiv <= ipwr[idec]) break;
  }
  jdec = 0; ibgn = idec + 1; ndbcd = (int) idpar[0]; nflg = -ndbcd;
  idcmp = &idpar[ndbcd * npar];
  for (k = 0; k < 3; ++k)
    idcmp[k + ndbcd] = (bcd_t) 0;
  /*      while there are unlabeled parameters */
  while (kbgn <= kend) {
    /*   open file of parameter names */
    if (jdiv <= 0) {
      lu = file_helpers::open_input_optional(pname);
      if (lu == NULL)
        break;
      /* first line contains divisor to ignore lower digits of ID */
      if (fgetstr(buf, NCARD, lu) <= 0)
        break;
      jdiv = atoi(buf);
      if (jdiv <= 0)
        jdiv = 1;
      for (jdec = 0; jdec < 3; ++jdec) {
        if (jdiv <= ipwr[jdec]) break;
      }
      jdec -= idec;
    }
    if (lu == NULL || fgetstr(buf, NCARD, lu) <= 0)
      break;
    /* parse into number + string */
    k = getbcd(buf, idcmp, nflg);
    while (buf[k] == ' ')  /* ignore spaces */
      ++k;
    parx = &buf[k];
    if (NEGBCD(idcmp[0]) != 0) {
      /*     negative ID signals that parameter name field is a new file */
      /*          to find further parameter names */
      fclose(lu);
      pname = parx;
      jdiv = 0;
    } else {
      if (k >= 82) /* check if string too short */
        continue;
      ibcd = kbgn * ndbcd;
      for (k = kbgn; k <= kend; ++k, ibcd += ndbcd) {
        for (i = ibgn; i < ndbcd; ++i) {
          if (idcmp[i + jdec] != idpar[i + ibcd]) break;
        }
        if (i != ndbcd) continue;
        pstr = &parlbl[k * nlbl];
        if (pstr[0] == '\0') {
          memcpy(pstr, parx, nlbl);
          if (k == kend) {
            --kend;
          } else if (k == kbgn) {
            ++kbgn;
          }
        }
      }
    }
  }
  if (lu != NULL)
    fclose(lu);
  return 0;
}                               /* getlbl */

int filbak(char *flu, char *fbak)
{
  /*     FLU= input file name (.par) */
  /*     FBAK =backup file name (.bak) */
  FILE *lu, *lubak;
  char cline[NCARD];

  lu = fopen(flu, "r");
  if (!lu) return 1;
  lubak = fopen(fbak, "w");
  if (!lubak) { fclose(lu); return 1; }
  while (fgets(cline, NCARD, lu)) {
    fputs(cline, lubak);
  }
  fclose(lu);
  rewind(lubak);  /* TODO: pointless? */
  fclose(lubak);
  return 0;
}                               /* filbak */

int prcorr(FILE *lufit, int nfit, double *cor, int ndcor, double *err, int normalize_to_correlation)
{ /*  calculates and prints correlation/covariance coefficients */
  int i, j, n, ndiag;
  double val, *dcor, *pcor;
  // Temporary array to store diagonal elements sqrt(Cov_ii) if normalizing
  double *sqrt_diag_cov = NULL; // Initialize to NULL
  double tiny = 10e-13;
  ndiag = ndcor + 1;
  dcor = cor;
  for (j = 0; j < nfit; ++j)
  { /* scale */
    if (j > 0)
      dcor += ndiag;
    n = nfit - j;
    val = 1. / err[j];
    dscal(n, val, dcor, 1);
  }
  for (n = 1; n < nfit; ++n)
  { /* matrix multiply */
    i = nfit - n;
    pcor = &cor[i];
    for (j = 0; j < i; ++j)
    {
      dcor[j - i] = ddot(n, pcor, 1, dcor, 1);
      pcor += ndcor;
    }
    dcor -= ndiag;
  }
  pcor = dcor = cor;
  for (n = 1; n < nfit; ++n)
  {
    pcor += ndcor;
    ++dcor;
    dcopy(n, pcor, 1, dcor, ndcor);
  }
  // --- Optional Normalization to True Correlation Matrix ---
  if (normalize_to_correlation && nfit > 0)
  {
    sqrt_diag_cov = (double *)malloc(nfit * sizeof(double));
    if (sqrt_diag_cov == NULL)
    {
      fprintf(lufit, "Error: Could not allocate memory for sqrt_diag_cov in prcorr.\n");
      // Proceed to print unnormalized covariance
    }
    else
    {
      // Calculate sqrt(Cov_ii)
      pcor = cor; // Points to cor[0,0]
      for (i = 0; i < nfit; ++i)
      {
        if (*pcor > 0.0)
        { // Cov_ii must be positive
          sqrt_diag_cov[i] = sqrt(*pcor);
        }
        else
        {
          sqrt_diag_cov[i] = 1.0; // Avoid division by zero, correlation will be 0 or undefined
        }
        pcor += ndiag; // Move to next diagonal Cov_i+1,i+1
      }

      // Normalize: Corr_ij = Cov_ij / (sqrt_diag_cov[i] * sqrt_diag_cov[j])
      dcor = cor; // To iterate through Cov matrix elements
      for (i = 0; i < nfit; ++i)
      {              // Row i
        pcor = dcor; // Start of row i (or column i)
        for (j = 0; j < nfit; ++j)
        { // Column j
          if (sqrt_diag_cov[i] != 0.0 && sqrt_diag_cov[j] != 0.0)
          {
            *pcor /= (sqrt_diag_cov[i] * sqrt_diag_cov[j]);
          }
          else
          {
            *pcor = 0.0; // Or NaN, if diagonal was zero
          }
          pcor++; // Next element in the row/column (assuming flat access for now)
                  // This needs to correctly access cor[i][j]
                  // If cor is col-major: cor[i + j*ndcor]
        }
        dcor += ndcor; // Next row/column start for outer loop
      }
      // More careful access for normalization:
      for (i = 0; i < nfit; ++i)
      { // Target row index
        for (j = 0; j < nfit; ++j)
        {                                                    // Target col index
          double *current_element_ptr = cor + j * ndcor + i; // Access cor[i,j] if col-major
          if (sqrt_diag_cov[i] > tiny && sqrt_diag_cov[j] > tiny)
          { // Avoid div by zero/small
            *current_element_ptr /= (sqrt_diag_cov[i] * sqrt_diag_cov[j]);
          }
          else
          {
            if (i == j)
              *current_element_ptr = 1.0; // Corr_ii = 1
            else
              *current_element_ptr = 0.0; // Undefined correlation
          }
        }
      }
    }
  }
  n = 0;
  dcor = cor;
  for (i = 1; i <= nfit; ++i)
  {
    pcor = dcor;
    for (j = 1; j <= nfit; ++j)
    {
      if (i != j)
      {
        fprintf(lufit, "%3d%3d%10.6f", i, j, *pcor);
        n = (n + 1) & 7;
        if (n == 0)
          fputc('\n', lufit);
      }
      ++pcor;
    }
    dcor += ndcor;
  }
  if (n > 0)
    fputc('\n', lufit);
  return 0;
} /* prcorr */

namespace {
  std::vector<SXLINE> g_lines;     // g_lines[0] = sentinel; g_lines[i] = line i (1-indexed)
  std::vector<double> g_dnudp_buf; // flat derivatives storage, size nlines * ndbl
}

SXLINE *line_at(int ipos) {
  return &g_lines[(size_t)abs(ipos)];
}

void init_line_buffer(size_t nlines, int nfit) {
  int ndbl = nfit > 0 ? nfit : 1;
  g_lines.assign(nlines + 1, SXLINE{});
  g_dnudp_buf.assign(nlines * (size_t)ndbl, 0.0);
  for (size_t i = 1; i <= nlines; ++i)
    g_lines[i].dnudp = &g_dnudp_buf[(i - 1) * (size_t)ndbl];
  g_lines[0].linku = g_lines[0].linkl = 0;
  g_lines[0].dnudp = nullptr;
}

void release_line_buffer() {
  g_lines.clear(); g_lines.shrink_to_fit();
  g_dnudp_buf.clear(); g_dnudp_buf.shrink_to_fit();
}

void dnuadd(int npar, int nparx, int initl, int indx, int ifac,
            double *egy, double *egyder, int nsize, int line,
            const double *par, const double *fac)
{                               /*     subroutine to add energies and derivatives to frequency lists */
  /*  ON INPUT: */
  /*     NPAR   = number of parameters */
  /*     NPARX  = number of parameters if INDX < 0 */
  /*     INITL  = first energy computed for this frequency  IF < 0 */
  /*     INDX   = index in energy list */
  /*     IFAC   = scaling for energies */
  /*     EGY    = energy vector */
  /*     EGYDER = energy derivative */
  /*     NSIZE   = column size of EGYDER */
  /*     LINE   = line counter for frequency derivatives */
  /*     PAR    = parameters */
  /*     FAC    = scaling values */
  SXLINE *xline;
  double f, dtmp, *deriv;
  long itmp;
  int noff, nparn;

  xline = line_at(line);
  deriv = xline->dnudp;
  if (initl < 0)
    xline->cfrq = 0.;
  f = fac[ifac];
  /* add derivatives */
  if (indx < 0) {               /* set up to ignore last few derivatives */
    indx = -1 - indx;
    nparn = nparx;
    noff = npar - nparn;
    itmp = nparn;
    itmp = indx + itmp * nsize;
    /*  add energies subtracting effect of ignored derivatives */
    dtmp = egy[indx] - ddot(noff, &egyder[itmp], nsize, par, 1);
    if (initl < 0)
      memset(&deriv[nparn], 0, noff * sizeof(double));
  } else {                      /* add energies */
    indx = indx - 1;
    dtmp = egy[indx];
    nparn = npar;
  }
  if (ifac != 0)
    dtmp *= f;
  xline->cfrq += dtmp;
  /*  sum derivatives */
  if (initl < 0) {
    dcopy(nparn, &egyder[indx], nsize, deriv, 1);
    if (ifac != 0)
      dscal(nparn, f, deriv, 1);
  } else {
    daxpy(nparn, f, &egyder[indx], nsize, deriv, 1);
  }
}                               /* dnuadd */

double dnuget(int iflg, int npar, double f, int line, double *dvec)
{
  SXLINE *xline;
  double frq;

  /*     subroutine to get derivatives and add to DVEC */
  /*  ON INPUT: */
  /*     IFLG = flag to determine type of operation          */
  /*     NPAR = number of parameters                         */
  /*     F    = scale factor (if IFLG >= 0 )                 */
  /*     LINE   = line counter for frequency derivatives     */
  /*  RETURNS: */
  /*     calculated frequency                                */
  /*     DVEC = modified vector of derivatives and frequency */

  xline = line_at(line);
  frq = xline->cfrq;
  if (iflg == 0) {
    /* copy scaled derivatives into DVEC */
    dcopy(npar, xline->dnudp, 1, dvec, 1);
    dscal(npar, f, dvec, 1);
    dvec[npar] = f * (xline->xfrq - frq);
  } else if (iflg > 0) {
    /*  sum scaled derivatives with DVEC */
    daxpy(npar, f, xline->dnudp, 1, dvec, 1);
    dvec[npar] += f * (xline->xfrq - frq);
  }
  return frq;
}                               /* dnuget */

int getdbk(int *link, int *iblk, int *indx, int *initl, int *ifac)
{
  SXLINE *xline;
  int iret, init, iup, ilow;

  /*   find data for line = abs(LINK)      neg sign points to lower */
  /*   returns II  = line number */
  /*           IBLK= requested block number */
  /*           INDX= requested index */
  /*           INITL < 0 if first block */
  /*           IFAC= code for energy conversion */
  /*                 0 for 1., 1 for -1., 2 for c, 3 for -c */
  iret = *link;
  xline = line_at(iret);
  iup = xline->ibu;
  ilow = xline->ibl;
  init = iup - ilow;
  if (iret >= 0) {
    if (ilow == 0 || init == 0)
      init = -1;
    *iblk = iup;
    *link = xline->linku;
    *indx = xline->inu;
    *ifac = 0;
  } else {
    iret = -iret;
    init = -init;
    *iblk = ilow;
    *link = xline->linkl;
    *indx = xline->inl;
    *ifac = 1;
  }
  *initl = init;
  /*  if XERR < 0 use cm-1 */
  if (xline->xerr < 0.)
    *ifac += 2;
  return iret;
}                               /* getdbk */

int frqdat(int line, int *ibln, double *txfrq, double *txwt, double *txerr, short *iqn)
{  /*  gets blend code, frequency, weight, error , and quantum numbers */
  static size_t ndqn = 2 * MAXQN * sizeof(short);
  SXLINE *xline;
  xline = line_at(line);
  *txfrq = xline->xfrq;
  *txerr = xline->xerr;
  *ibln = xline->bln;
  *txwt = xline->xwt;
  memcpy(iqn, xline->qn, ndqn);
  return xline->ibu;
}                               /* frqdat */

int lnlink(int *prvblk, int nblk, int iblk, int line)
{
  int i, last;
  SXLINE *xline;
  /*   find place to insert link for line */
  /*   return pointer to next line */
  /*   PRVBLK contains pointer to last line using given block */

  if (iblk <= 0)
    return 0;
  if (iblk > nblk)
    iblk = nblk;
  last = 0;
  for (i = iblk; i >= 0; --i) {
   /* search down for non zero element of PRVBLK */
    last = prvblk[i];
    if (last != 0)
      break;
  }
  prvblk[iblk] = line;
  xline = line_at(last);
  if (last >= 0) {
    last = xline->linku;
    xline->linku = line;
  } else {
    last = xline->linkl;
    xline->linkl = line;
  }
  if (line != 0) {
    xline = line_at(line);
    if (line >= 0) {
      xline->linku = last;
    } else {
      xline->linkl = last;
    }
  }
  return last;
}                               /* lnlink */
