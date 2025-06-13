/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   25 March 1999: read option cards with fgetstr */
/*   30 Dec.  1999: include changes for dlsq */
/*   10 Oct.  2001: change fit diverging code */
/*   21 Sept. 2002: fix NRJ criterion */
/*   18 Aug.  2003: code cleanup, @comment@ is for splint */

/**************************************************************************/
/*                                                                        */
/*   THIS IS A GENERALIZED LINE FITTING PROGRAM                           */
/*   IT FITS LINES TO PARAMETERS IN A MODEL HAMILTONIAN                   */
/*   BLENDED LINES ARE TREATED SPECIALLY:                                 */
/*     IF THE EXPTL.FREQ. IS THE SAME TO 1 HZ THEN ONLY THE INVERSE ERROR */
/*             AVERAGED FREQUENCIES ARE USED IN THE FIT FOR THE BLEND     */
/*                                                                        */
/**************************************************************************/
/* "trust region" Marquardt fitting is described in John. E. Dennis and   */
/* Robert B. Schnabel, Numerical Methods for Unsconstrained Optimization  */
/* and Non-linear Equations, Prentice-Hall, 1983.                         */

#include <stdio.h>
#include <string.h>
#include <string>
#include <math.h>
#include <stdlib.h>
//#include "calpgm.h"
#include "lsqfit.h"
#include "SpinvEngine.hpp"
#include "DpiEngine.hpp"
#include "CalFit.hpp"

#define NDCARD 130
#define PR_DELAY 6    /* seconds delay between informational messages */
/************** CALFIT interfaces ***********************************/
/**
 * @brief Format quantum numbers for output
 *
 * @param nqn Number of quantum numbers
 * @param qnum Array of quantum numbers
 * @param aqnum Output string for formatted quantum numbers
 * @return int Always returns 0
 */
int qnfmt2(int nqn, short *qnum, /*@out@*/ char *aqnum);

/**
 * @brief Format parameter values, errors, and changes for output
 *
 * @param par Parameter value
 * @param errx Parameter error
 * @param dif Parameter change
 * @param ptmp Output string for formatted parameter
 * @return int Always returns 0
 */
int parer(double par, double errx, double dif,
                 /*@out@*/ char *ptmp);

/**
 * @brief Read experimental lines from input file
 *
 * @param luin Input file pointer
 * @param nline Pointer to number of lines
 * @param iqnfmt Quantum number format for line input
 * @return int Largest quantum number encountered
 */
int linein(FILE * luin, int *nline, int iqnfmt);

/**
 * @brief Process lines and set up block structure for fitting
 *
 * @param calc Engine used for finding relevant quantum numbers
 * @param lu Output file pointer for listing
 * @param flg Flag for printing (negative for detailed output)
 * @param nline Number of lines
 * @param nblkpf Number of blocks per F quantum number
 * @param iqnfmt Quantum number format for line input
 * @return int Number of bad lines
 */
int lineix(CalculationEngine *calc, FILE *lu, int flg,
           int nline, int nblkpf, int iqnfmt);

int getblk(CalculationEngine *calc, /*@out@*/ int *iblk,
           /*@out@*/ int *indx, short *iqnum,
           int nblkpf, int ipos, int nqn);
/********************************************************************/
static char card[NDCARD];  /* Buffer for reading input lines */

/**
 * @brief Main function for SPFIT - Spectroscopic Parameter Fitting program
 *
 * This program fits spectroscopic parameters to experimental line frequencies.
 * It uses the CalFit class to perform the fitting process, delegating the core
 * logic to a modular implementation.
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line arguments
 * @return int Exit status (0 for success)
 */
int main(int argc, char *argv[])
{
#define NFILE 5
  static const char *ext[NFILE] = { "par", "lin", "fit", "bak", "var" };
  enum efile {epar, elin, efit, ebak, evar};
  char *fname[NFILE+1];
  std::string engineType = "spinv";

  // Choose the calculation engine based on the first argument
  if (argc > 1) {
    if (strcasecmp(argv[1], "--dpi") == 0) {
      engineType = "dpi";
      char *name = *argv; argv++; *argv = name; argc--; // Remove first argument
    } else if (strcasecmp(argv[1], "--spinv") == 0) {
      char *name = *argv; argv++; *argv = name; argc--; // Remove first argument
    }
  }

  // Open read and write files
  filget(argc, argv, NFILE, fname, ext);

  // Create CalFit instance with the selected engine
  CalFit calFit(engineType);

  // Prepare input and output structures
  CalFitInput input;
  CalFitOutput output;

  // Read input data using CalFitIO
  if (!CalFitIO::readInput(fname[epar], fname[elin], input)) {
    puts("Failed to read input files");
    return EXIT_FAILURE;
  }

  // Run the fitting process
  if (!calFit.run(input, output)) {
    puts("Fitting process failed");
    return EXIT_FAILURE;
  }

  // Write output data using CalFitIO
  if (!CalFitIO::writeOutput(fname[efit], fname[ebak], fname[evar], output)) {
    puts("Failed to write output files");
    return EXIT_FAILURE;
  }

  puts("FIT COMPLETE");
  return 0;
}                               /* MAIN */

/**
 * @brief Format quantum numbers as a string for output
 *
 * Converts an array of quantum numbers to a formatted string representation.
 * Each quantum number is formatted as a 3-character field, and the string
 * is padded to accommodate up to 12 quantum numbers.
 *
 * @param nqn Number of quantum numbers to format
 * @param qnum Array of quantum numbers
 * @param aqnum Output string buffer for formatted quantum numbers
 * @return int Always returns 0
 */
int qnfmt2(int nqn, short *qnum, char *aqnum)
{
  /* Local variables */
  int i;

  /*  formats quantum numbers for output */
  /*     NQN   = number of quanta */
  /*     QNUM  = vector of quanta */
  /*     AQNUM = string of quanta in character form */
  for (i = 0; i < nqn; ++i) {
    sprintf(aqnum, "%3d", (int) qnum[i]);
    aqnum += 3;
  }
  for (i = nqn; i < 12; ++i) {
    aqnum[2] = aqnum[1] = aqnum[0] = ' ';
    aqnum += 3;
  }
  aqnum[0] = '\0';
  return 0;
}                               /* qnfmt2 */

/**
 * @brief Format parameter values, errors, and changes for output
 *
 * Creates a specially formatted string representation of a parameter value,
 * its error, and its change. The function automatically determines the
 * appropriate number format (fixed or scientific notation) based on the
 * magnitude of the values.
 *
 * @param par Parameter value
 * @param errx Parameter error
 * @param dif Parameter change
 * @param ptmp Output string buffer for formatted parameter
 * @return int Always returns 0
 */
int parer(double par, double errx, double dif, char *ptmp)
{
  static int czero = (int) '0';  /* ASCII code for '0' character */
  char *pfmt;                    /* Pointer for building format string */
  double adif, apar, aten, aerr; /* Working copies of input values */
  char chexp[6], fmt[34];        /* Format string components */
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
  ie = (int) (log10(fabs(aerr) + 1.e-37) - 102.5) + 100;
  id = (int) (log10(fabs(adif) + 1.e-37) - 100.0) + 100;
  ip = (int) (log10(fabs(apar) + 1.e-37) - 100.0) + 100;
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
  if (msd <= -2) {              /* number too small without exponent */
    k = (1 - msd) / 3;
    efield = -3 * k;
    while ((--k) >= 0)
      aten *= 1000;
  } else if (lsd < 0) {         /* number too big without exponent */
    k = (1 + msd) / 3;
    if (k > 0)
      efield = 3 * k;
    while ((--k) >= 0)
      aten *= 0.001;
  }
  if (efield != 0) {            /* E format */
    lsd += efield;
    memcpy(chexp, "0fE+00", 6);
    if (efield < 0) {
      chexp[3] = '-';
      efield = -efield;
    }
    msd = efield / 10;
    if (msd > 0) {
      efield -= msd * 10;
      chexp[4] = (char) (msd + czero);
    }
    chexp[5] = (char) (efield + czero);
    apar *= aten;
    aerr *= aten;
    adif *= aten;
  } else {                      /* F format */
    memcpy(chexp, "0f    ", 6);
  }
  if (lsd > 9)
    lsd = 9;
  if (lsd > 0)
    chexp[0] = (char) (lsd + czero);
  while ((lsd--) > 0)
    aerr *= 10.;
  ie = (int) (aerr + 0.5);
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
}                               /* parer */


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
 * @param iqnfmt Quantum number format for line input
 * @return int Largest quantum number encountered
 */
int linein(FILE *luin, int *nline, int iqnfmt)
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
  for (i = 1; i <= mxline; ++i) {       /*  loop for reading lines */
    xline = lbufof(1, i);
    iqnum = xline->qn;
    if (getlin(luin, nqn, nqnt, iqnum, &xfrqn, &xerrn, &xwtn,
               card, NDCARD) < 0) {
      *nline = i - 1;
      return mxqn;
    }
    iqf = iqnum[nqnu];
    if (iqf == -1) {
      if (kqnu >= 0) {
        iqf = -iqnum[kqnu];
        if (iqf >= 0)
          iqf = -1;
      }
      iqnum[nqnu] = (short) iqf;
    }
    if (iqf < 0)
      iqf = -iqf;
    if (iqf > mxqn)
      mxqn = iqf;
    iqf = iqnum[nqnl];
    if (iqf == -1) {
      if (kqnl > 0) {
        iqf = -iqnum[kqnl];
        if (iqf >= 0)
          iqf = -1;
      }
      iqnum[nqnl] = (short) iqf;
    }
    if (iqf < 0)
      iqf = -iqf;
    if (mxqn < iqf)
      mxqn = iqf;
    xline->xfrq = xfrqn;
    xline->xerr = (float) xerrn;
    xline->xwt = (float) fabs(xwtn);
    xline->linku = 0;
    xline->linkl = 0;
    isblnd = 0;
    if (icmp != 0 && fabs(xfrqn - xfrqx) < fabs(xfrqn) * 1.e-14 + 1.e-8) {
      /* frq match */
      if (fabs(xerrn - xerrx) < 1e-7) {
        isblnd = 1;
      } else if ((xerrn / xerrx) > 2.0 && nbln > 2) {
        isblnd = 1; ++nbln; icmp = 0;
        xline->xwt = (float)0.;
        iqnum[0] = (short)-1;
        iqnum[nqn] = iqnum[0];
      }
    }
    if (isblnd != 0) {
      xline->bln = nbln;
      xline = lbufof(1, i - 1);
      xline->bln = -2;
      nbln += 2;
    } else {
      xline->bln = 0;
      nbln = 2; icmp = 1;
    }
    if (ipace <= i) {
      ipace += 100;
      printf("Reading Line %d\n", i);
      fflush(stdout);
    }
    xerrx = xerrn;
    xfrqx = xfrqn;
  }
  return mxqn;
}                               /* linein */

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
int lineix(CalculationEngine *calc, FILE *lu, int flg, int nline, int nblkpf, int iqnfmt)
{  /*   get lines from input and store them */
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
  char aqnum[6*MAXQN+2];

  nbad = 0;
  nblk = 0;
  prvblk = (int *) mallocq((size_t) (nsort + 1) * sizeof(int));
  prvblk[0] = 0;
  for (i = 1; i <= nsort; ++i) {
    prvblk[i] = 0;
  }
  nqn = iqnfmt % 10;
  if (nqn == 0) nqn = 10;
  nqn2 = nqn + nqn; ncat = nqn2;
  if (ncat < 12) ncat = 12;
  i = (iqnfmt / 100) % 5;
  if (i >= nqn) {
    ipos = 1;
  } else {
    ipos = nqn;
  }
  if (flg < 0) {
    fputs(" LINE,BLKU,INDXU,BLKL,INDXL,QUANTUM NUMBERS", lu);
    for (i = 0; i < 19; ++i)
      fputc(' ', lu);
    fputs("ENERGY    EXP. ERROR    WEIGHTS\n", lu);
  }
  xnorm = 0.;
  ipace = 50;
  /*       loop for converting lines */
  for (nread = 1; nread <= nline; ++nread) {
    xline = lbufof(1, nread);
    xfrqn = xline->xfrq;
    xerrn = xline->xerr;
    xwtn = xline->xwt;
    /* find blocks and index for upper and lower states */
    iqnum = xline->qn;
    getblk(calc, &iblku, &indxu, iqnum, nblkpf, ipos, nqn);
    getblk(calc, &iblkl, &indxl, &iqnum[nqn], nblkpf, ipos, nqn);
    if (iblkl == 0 && iqnum[nqn] >= 0)
      iblku = 0;
    xline->ibu = iblku;
    xline->inu = (short) indxu;
    xline->ibl = iblkl;
    xline->inl = (short) indxl;
    if (iblku == 0 && (xline->bln & 1) == 0) {
      /*  print out bad line and try for next */
      ++nbad;
      xline->xwt = 0.;
      xwtn = 0.;
      qnfmt2(nqn2, iqnum, aqnum);
      printf(    "Bad Line(%3d): %s %14.5f %8.5f\n",
                 nread, aqnum, xfrqn, xerrn);
      fprintf(lu,"Bad Line(%3d): %s %14.5f %8.5f\n",
                 nread, aqnum, xfrqn, xerrn);
    } else {
      /*  set up links for calculating in order of block */
      if (iblku <= iblkl) {
        lnlink(prvblk, nsort, iblku, nread);
        lnlink(prvblk, nsort, iblkl, -nread);
        if (nblk < iblkl)
          nblk = iblkl;
      } else {
        lnlink(prvblk, nsort, iblkl, -nread);
        lnlink(prvblk, nsort, iblku, nread);
        if (nblk < iblku)
          nblk = iblku;
      }
      if (flg < 0) {
        iqnum = xline->qn;
        fprintf(lu," %4d%4d%4d%4d%4d:", nread, iblku, indxu, iblkl, indxl);
        for (i = 0; i < ncat; ++i) {
          j = iqnum[i];
          fprintf(lu, "%3d", j);
        }
        fprintf(lu, " %14.4f %9.4f %9.4f", xfrqn, xerrn, xwtn);
        j = xline->bln;
        if (j != 0) {
          fprintf(lu, "   Line Blended with %3d\n", nread - (j >> 1));
        } else {
          fputc('\n', lu);
        }
      }
    }
    /* let the user know something is happening */
    if (ipace <= nread || nread == nline) {
      ipace += 50;
      printf("Converting Line %d\n", nread);
      fflush(stdout);
    }
    j = xline->bln;
    if (j != 0) {
      xnorm += xwtn;
      if (j > 0) {     /* normalize weights */
        xnorm = 1. / xnorm;
        for (j = nread - (j >> 1); j <= nread; ++j) {
          xline = lbufof(1, j);
          xline->xwt = (float) (xline->xwt * xnorm);
        }
        xnorm = 0.;
      }
    } else {
      xline->xwt = 1.;
    }
  }                             /* end loop for converting lines */
  orgblk = nsort;
  while (nblk > orgblk) {       /* finish up links */
    --orgblk;
    linkx = 0;
    for (j = 0; j < nsort; ++j) {
      if (prvblk[j] != 0) {
        linkx = prvblk[j];
        prvblk[j] = 0;
      }
    }
    prvblk[0] = linkx;
    prvblk[nsort] = 0;
    linkx = lnlink(prvblk, nsort, 1, 0);
    while (linkx != 0) {
      linky = linkx;
      getdbk(&linkx, &iblkl, &j, &j, &j);
      iblkl -= orgblk;
      lnlink(prvblk, nsort, iblkl, linky);
    }
    orgblk += nsort;
  }
  free(prvblk);
  return nbad;
}                               /* lineix */

int getblk(CalculationEngine *calc, int *iblk, int *indx,
           short *iqnum, int nblkpf, int ipos, int nqn)
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
