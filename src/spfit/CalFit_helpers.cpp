/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include <cstring> // for memset
#include <cmath>   // for fabs, sqrt
#include "splib/lsqfit.h" // for linear algebra functions and other definitions
#include <stdio.h>  // for FILE and fprintf
#include <string.h> // for strcpy
#include <vector>   // for std::vector
#include <algorithm>// for std::copy
#include "CalFitIO.hpp"
#include "CalFit.hpp"
#include "splib/lsqfit.h"
#include "splib/calpgm_types.h"
#include "splib/ulib.h"
#include "spfit/subfit.h"
#include "common/CalError.hpp"
#include "engine/SpinvEngine.hpp"
#include "engine/DpiEngine.hpp"

// original static functions now integrated as private members of class

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
  {
    // Pad to 12 quantum numbers (36 chars)
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
 * @param ptmp Output string for formatted parameter
 * @return int Always returns 0
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
 *
 * Reads spectral line data from the input file, including frequencies,
 * errors, weights, and quantum numbers. Handles blended lines by
 * identifying lines with matching frequencies. Stores the data in
 * the global line buffer for later processing.
 *
 * @param luin Input file pointer
 * @param nline Pointer to number of lines (input: max lines, output: actual lines)
 * @param iqnfmt Quantum number format for line input (from setfmt)
 * @return intLargest quantum number encountered
 */
/**
 * Parse a pre-read line card (equivalent to getlin without fgetstr).
 * Returns nc on success, -1 on parse failure, -2 on ncard too small.
 */
static int getlin_str(const char *line, int line_len, int nqn, short *idqn,
                      short *iqn, double *xfreq, double *xerr, double *xwt,
                      char *card, int ncard)
{
  static double terr = 1.e-07;
  static double txwt = 1.e-30;
  static double big = 1000.;
  static int fmt[20] = {3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3};
  double val[20];
  int j, jd, k, nn, nc;
  char ctmp;

  if (nqn <= 6) { nn = 12; nc = 36; }
  else { nn = nqn + nqn; nc = nn * 3; }
  if (ncard <= nc) return -2;
  if (line_len < nc) return -1;
  strncpy(card, line, (size_t)(ncard - 1));
  card[ncard - 1] = '\0';

  dcopy(nn, &big, 0, val, 1);
  ctmp = card[nc]; card[nc] = '\0';
  pcard(card, val, nn, fmt);
  card[nc] = ctmp;
  for (k = 0; k < nn; ++k) {
    j = (int)val[k];
    if (j > 999) { j = 0; jd = idqn[k]; if (jd >= 0) j = iqn[jd]; }
    iqn[k] = (short)j;
  }
  val[0] = 0.; val[1] = 0.01; val[2] = 1.;
  if (pcard(card + nc, val, 3, NULL) == 0) return -1;
  *xfreq = val[0];
  *xerr = val[1];
  if (fabs(*xerr) < terr) *xerr = terr;
  *xwt = val[2];
  if (*xwt < txwt) *xwt = txwt;
  return nc;
}

int CalFit::linein(const std::vector<std::string> &lines, int *nline, int iqnfmt)
{
  /* Local variables */
  SXLINE *xline;
  double xfrqn, xerrn, xwtn;
  double xfrqx, xerrx;
  int nqn, nqnu, nqnl;
  int kqnu, kqnl;
  int i, iqf, ipace, mxline;
  int mxqn, isblnd, icmp;
  short nbln, nqnt[20], *iqnum;
  char card[NDCARD];

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
  size_t line_idx = 0;
  for (i = 1; i <= mxline; ++i)
  { /*  loop for reading lines */
    xline = line_at(i);
    iqnum = xline->qn;
    if (line_idx >= lines.size() ||
        getlin_str(lines[line_idx].c_str(), (int)lines[line_idx].size(),
                   nqn, nqnt, iqnum, &xfrqn, &xerrn, &xwtn,
                   card, NDCARD) < 0)
    {
      *nline = i - 1;
      return mxqn;
    }
    ++line_idx;
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
      if (fabs(xerrn - xerrx) < 1e-7) {
        isblnd = 1;
      }
      else if ((xerrn / xerrx) > 2.0 && nbln > 2) { // TODO: should ensure xerrx is not zero
        isblnd = 1;
        ++nbln;
        icmp = 0;
        xline->xwt = (float)0.;
        iqnum[0] = (short)-1;
        iqnum[nqn] = iqnum[0];  // TODO: check if nqn>0 ?
      }
    }
    if (isblnd != 0)
    {
      xline->bln = nbln;
      xline = line_at(i - 1);
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
 *
 * Processes the spectral lines read by linein(), determining the Hamiltonian
 * blocks and indices for upper and lower states of each transition.
 * Sets up links for calculating energies and derivatives in block order.
 * Identifies and reports bad lines (those with invalid quantum numbers).
 *
 * @param lu Output file pointer for listing
 * @param flg Flag for printing (negative for detailed output)
 * @param nline Number of lines
 * @param nblkpf Number of blocks per F quantum number
 * @param iqnfmt Quantum number format for line input
 * @return int Number of bad lines
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
  int nblk, ipos_for_getblk, i, j, ipace, nread, iblkl, iblku, ncat;
  int linkx, indxl, linky, indxu, orgblk, nqn_for_hamiltonian_basis, nqn2, nbad;
  /*@owned@*/ int *prvblk;
  short *iqnum;
  char aqnum[6 * MAXQN + 2];

  nbad = 0;
  nblk = 0;
  prvblk = (int *)calalloc((size_t)(nsort + 1) * sizeof(int));
  prvblk[0] = 0;
  for (i = 1; i <= nsort; ++i)
  {
    prvblk[i] = 0;
  }
  nqn_for_hamiltonian_basis = iqnfmt % 10;
  if (nqn_for_hamiltonian_basis == 0)
    nqn_for_hamiltonian_basis = 10;
  nqn2 = nqn_for_hamiltonian_basis + nqn_for_hamiltonian_basis;
  ncat = nqn2;
  if (ncat < 12)
    ncat = 12;
  i = (iqnfmt / 100) % 5;
  if (i >= nqn_for_hamiltonian_basis)
  {
    ipos_for_getblk = 1;
  }
  else
  {
    ipos_for_getblk = nqn_for_hamiltonian_basis;
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
    xline = line_at(nread);
    xfrqn = xline->xfrq;
    xerrn = xline->xerr;
    xwtn = xline->xwt;
    /* find blocks and index for upper and lower states */
    iqnum = xline->qn;
    this->getblk(&iblku, &indxu, iqnum, nblkpf, ipos_for_getblk, nqn_for_hamiltonian_basis);
    this->getblk(&iblkl, &indxl, &iqnum[nqn_for_hamiltonian_basis], nblkpf, ipos_for_getblk, nqn_for_hamiltonian_basis);
    if (iblkl == 0 && iqnum[nqn_for_hamiltonian_basis] >= 0)
    {
      iblku = 0;
    }
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
      {
        /* normalize weights */
        xnorm = 1. / xnorm;
        for (j = nread - (j >> 1); j <= nread; ++j)
        {
          xline = line_at(j);
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
  {
    /* finish up links */
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
} /* getblk */
