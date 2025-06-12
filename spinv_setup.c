/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   24 March  1999: Fixed sin scaling in dircos */
/*   25 March  1999: Fixed diag = 2 error in call to ordham*/
/*   27 March  1999: Changed diag = 3 to take care of l-doubled states */
/*   30 March  1999: Fixed sin scaling bug */
/*   21 April  1999: Fixed bug at end of hamx for KROLL case */
/*   22 April  1999: Fixed bug in Itot option */
/*   25 May    1999: fix bug when largest number of spins is not default */
/*   23 July   1999: fix bug in SETOPT that returned incorrect IQNFMT */
/*    7 Sept.  1999: fix bug in GETQQ for that filled ivs with short */
/*   29 Dec.   1999: improve ORDHAM */
/*   12 April  2000: negate Hamiltonian for kappa > 0 on first pass */
/*   12 April  2000: separate ORDHAM and FIXHAM */
/*   12 April  2000: subtract constant from energy for more precision*/
/*    1 August 2000: special code for unit matrix */
/*   12 Sept.  2000: fix up phase for cross-multiplet interaction */
/*    9 Oct.   2000: change dipole type 10 to double commutator with N^2 */
/*   29 Nov.   2000: redefine idpars and spar to put more calc. in initial */
/*    5 Feb.   2001: define augmented Euler and Fourier series */
/*    7 July   2001: fix bug for l doublet pairs in setopt */
/*   26 Sept.  2001: fix bug in idpari for oblate 6x and 2x interactions */
/*   12 March  2002: improve selections for off-diagonal l in dircos */
/*   27 June   2002: increase vibration limit to 359 */
/*   29 July   2002: fix bug in fixham (DIAG = 2,3)*/
/*   18 Aug.   2003: code cleanup, @comment@ is for splint */
/*   24 Nov.   2003: implement code for higher symmetry IAX options */
/*    3 Feb    2004: implement code for higher symmetry spins */
/*   19 May    2004: implement code for bcd storage of idpar and idip */
/*   20 July   2004: fix bug in code for bcd storage of idip */
/*   19 Jan.   2005: extended Euler operator */
/*   28 May    2005: fix bug in getmask for operator between l-doubled state */
/*   27 Sept.  2005: fix bug for 6-fold symmetry B-weights in getqn & setopt */
/*   17 Oct.   2005: fix bug in getmask for operator within l-doubled state */
/*   12 Dec.   2005: fix bug Itot basis and include mixed basis */
/*   10 Sept.  2006: add other phase conventions */

/* split from original spinv.c */

/* Functions responsible for initializing the program state, reading options, and setting up primary data structures based on input files. */
/* Also includes definitions / initialization of all externed variables in spinv_internal.h */
/*
setopt(): Reads and processes the main option card(s) from the .par file.
setblk(): Initializes Hamiltonian block structures after options and parameters are read.
pasort(): Parses and sorts parameter IDs, populating the spar_head linked lists.
setfmt(): Determines quantum number output formats.
setint(): Initializes dipole data from the .int file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calpgm.h"
#include "spinit.h"
#include "spinv_internal.h"
#include "SpinvContext.hpp"


/**
 * @brief Initializes dipole data from BCD encoded dipole identifiers.
 *        Parses each dipole ID, determines its type, symmetry, involved vibrational states,
 *        and spin operators. It also sets flags for imaginary dipoles based on operator
 *        order and selected phase conventions. Invalid or unsupported dipole types are flagged.
 * @param list_unit File pointer for writing warning messages.
 * @param ifdiag_out Output: Pointer to a BOOL that will be set to TRUE if the global diagonalization
 *                   option (glob.idiag) is non-negative, FALSE otherwise.
 * @param nsav_out Output: Pointer to an integer, always set to 1 in this function.
 * @param num_dipole_params Total number of dipole parameters to process.
 * @param bcd_dipole_ids Array of dipole identifiers in BCD format. The first element bcd_dipole_ids[0]
 *                       is expected to contain the number of BCD digits per ID.
 * @param is_imaginary_dipole_out Output array of integers, where each element is set to 1 if the
 *                                corresponding dipole is determined to be imaginary based on phase
 *                                conventions, 0 otherwise.
 * @return int Always returns 0 after processing.
 */
int setint(struct SpinvContext* ctx, FILE *lu, BOOL *ifdiag, int *nsav, const int ndip, bcd_t *idip, int *isimag)
{ /*   subroutine to initialize intensity data */
  static char *bad_msg[] = {"\n", "(1) bad symmetry field\n",
                            "(2) vib out of range\n", "(3) unknown  dipole type\n",
                            "(4) transition between states of different weight\n",
                            "(5) alpha is not 0\n", "(6) bad spin field\n",
                            "(7) K is greater tnan Kmax\n"};
  SVIB *pvib1, *pvib2;
  short *iiv1, *iiv2;
  double zfac;
  unsigned int ijv;
  int i, icase, isym, ivdif, ii, iv1, iv2, lv1, lv2, ifac, lv;
  int si1, ibcd, ndecv, nbcd, iflg, kd, ld, ifc, imag, ldel;
  int nimag[16], ioff, kmin, nmin, nn, kavg, good_case, bad_dip;
  bcd_t itmp;
  *ifdiag = (ctx->glob.idiag >= 0);
  *nsav = 1;
  if (ctx->nddip > 0)
  {
    free(ctx->dipinfo);
    ctx->dipinfo = NULL;
  }
  ctx->dipinfo = (SDIP *)mallocq((size_t)ndip * sizeof(SDIP));
  ctx->dipinfo0.fac = 0.;
  memcpy(ctx->dipinfo, &ctx->dipinfo0, sizeof(SDIP));
  ctx->nddip = ndip;
  for (i = 0; i < 16; ++i)
  {
    nimag[i] = 0;
  }
  ifac = ctx->glob.vibfac + 1;
  ndecv = ctx->glob.vibdec;
  nbcd = (int)idip[0] & 0x7f;
  for (i = 0, ibcd = 0; i < ndip; ++i, ibcd += nbcd)
  {
    ii = (int)idip[ibcd + ndecv + 1];
    si1 = (ii >> 4) & 0x0f;
    ii &= 0x0f;
    if (ndecv > 2)
      ii = 100 * ii + bcd2i(idip[ibcd + 3]);
    if (ndecv > 1)
      ii = 100 * ii + bcd2i(idip[ibcd + 2]);
    itmp = idip[ibcd + 1];
    isym = (int)itmp & 0x0f;
    ii = 10 * ii + ((int)(itmp >> 4) & 0x0f);
    iv2 = ii / ifac;
    iv1 = ii - iv2 * ifac;
    ii = (int)idip[ibcd + ndecv + 2];
    icase = ii & 0x0f;
    ifc = (ii >> 4) & 0x0f;
    if (icase == 7 || icase == 8)
    {
      ifc = 10 * ifc + si1;
      si1 = 0;
    }
    else
    {
      icase += ifc * 10;
      ifc = 0;
    }
    kavg = 0;
    ivdif = iv1 - iv2;
    if (ivdif < 0)
    {
      ii = iv1;
      iv1 = iv2;
      iv2 = ii;
    }
    pvib1 = &ctx->vinfo[iv1];
    pvib2 = &ctx->vinfo[iv2];
    lv1 = pvib1->lvupper;
    lv2 = pvib2->lvupper;
    if (ODD(lv1))
      --iv1;
    if (ODD(lv2))
      --iv2;
    lv1 = pvib1->lvqn;
    lv2 = pvib2->lvqn;
    /* MAKE SURE VIBRATIONS ARE ORDERED */
    iflg = 1;
    bad_dip = 0;
    zfac = 1.;
    do
    { /* one pass */
      /*  CONVERT SYMMETRY */
      bad_dip = 1;
      if (isym > 3 || icase > MAXINT)
        break;
      bad_dip = 2;
      if (iv1 >= ctx->glob.nvib)
        break;
      if (iv2 >= ctx->glob.nvib)
        break;
      if (ctx->glob.oblate == FALSE)
        isym = ctx->revsym[isym];
      if (isym == 0)
      {
        ld = 0;
        zfac = 0.009274;
        ioff = 8;
      }
      else
      {
        ld = 1;
        iflg |= MELEC;
        ioff = 0;
      }
      imag = 0;
      kd = ctx->isoddk[isym]; /* ld = 0,1; kd = 0, 1, 1, 0 */
      good_case = 4;
      if (icase > 0)
      {
        switch (icase)
        {
        case 1: /* ld = 2,2; kd = 0, 1, 1, 2 */
          ld = 2;
          if (isym == 3)
            kd = 2;
          break;
        case 2: /* ld = 2,2; kd = 2, 1, 1, 2 */
        case 11:
        case 12:
          ld = 2;
          if (kd == 0)
            kd = 2;
          if (icase == 11)
          {
            iflg |= MINOQ;
            ++imag;
            good_case = isym; /* commutator */
          }
          else if (icase == 12)
          {
            ifc = -1;
            ++imag;
            isym = 3 - isym;
          }
          break;
        case 5: /* ld = 1,2; kd = 0, 1,  1, 0 */
          if (isym == 0)
          {
            iflg |= MLZ;
            ld = 1;
            ++imag;
            if (lv1 != lv2 || lv1 <= 0 || si1 > 0)
              good_case = 0;
          }
          else
          {
            iflg |= MINOQ;
            ++imag;
            good_case = isym; /* commutator */
          }
          break;
        case 6: /* ld = 0,3; kd = 2, 3, 3, 2 */
          if (isym != 0)
            ld = 3;
          kd += 2;
          break;
        case 7:
          break;
        case 8:
          ifc = -1 - ifc;
          ++imag;
          isym = 3 - isym;
          break;
        case 10:
          iflg |= MINOQ;
          good_case = isym; /* double commutator */
          break;
        }
      }
      bad_dip = 3;
      if (good_case == 0)
        break;
      imag += ld;
      /* nsym = */ setgsym(ctx, (int)pvib1->gsym);
      bad_dip = 4;
      if (testwt(ctx, pvib1, pvib2, isym, 0))
        break;
      ldel = lv1 - lv2;
      if (ldel != 0)
      {
        bad_dip = 3;
        if (icase == 5 && isym == 0)
          break;
        if (ldel < 0)
          ldel = -ldel;
        if (ivdif < 0) /* DELTA K * DELTA L < 0 */
          ldel = -ldel;
      }
      else if (lv1 != 0)
      {
        if (ctx->glob.newlz)
        {
          if (lv1 < 0)
          {
            iflg |= MLZ;
            isym = 3 - isym;
            ++imag;
          }
        }
        else
        {
          if (lv1 < 0)
            break;
        }
      }
      if (ODD(imag))
        iflg |= MODD;
      /* check for imaginary operators */
      if (iv1 == iv2 && ldel == 0)
      {
        if (ODD(imag))
        {
          if (isym == 0)
            break;
          nimag[isym + ioff] += 1;
          iflg |= MIMAG;
        }
        else if (isym != 0)
        {
          nimag[isym + ioff + 4] += 1;
          iflg |= MIMAG;
        }
      }
      if (ctx->glob.nitot >= 3)
      {
        bad_dip = 5;
        if (MOD(kd - ldel, ctx->glob.nitot) != 0)
          break; /* require alpha == 0 */
        if (si1 <= ctx->itsym)
          iflg |= MIDEN; /* identity operator */
      }
      lv = 1;
      if (si1 > 0)
      {
        /* check spin properties */
        bad_dip = 6;
        if (ld == 0 && isym != 0)
          break;
        iiv1 = pvib1->spt;
        iiv2 = pvib2->spt;
        if (checksp(ctx, TRUE, si1, 0, iiv1, iiv2, &zfac) != 0)
          break;
        lv = TEST(iflg, MELEC) ? 2 : 0; /*  simple electric dipole */
      }
      if (ifc < 0)
        iflg |= MFCM;
      if (ctx->glob.nofc && ifc != 0)
      {
        if (ifc < 0)
          ifc = -1 - ifc;
        kavg = ifc;
        ifc = 0;
        bad_dip = 7;
        if (kavg > pvib1->knmax && kavg > pvib2->knmax)
          break;
      }
      if (ld == kd && ld < lv)
        ld -= (ld - lv) & 0x7fe; /* subtract even number */
      ctx->dipinfo[i].fac = zfac;
      ctx->dipinfo[i].flg = (short)iflg;
      ctx->dipinfo[i].kd = (signed char)kd;
      ctx->dipinfo[i].ld = (signed char)ld;
      ctx->dipinfo[i].ldel = (signed char)ldel;
      ctx->dipinfo[i].fc = (signed char)ifc;
      ctx->dipinfo[i].kavg = (char)kavg;
      bad_dip = 0;
    } while (FALSE);
    if (bad_dip != 0)
    {
      if (bad_dip > 0)
      {
        putbcd(ctx->sbcd, NSBCD, &idip[ibcd]);
        if (bad_dip * sizeof(char *) >= sizeof(bad_msg))
          bad_dip = 0;
        fprintf(lu,
                " WARNING: dipole %3d %s has no matrix elements, %s",
                (i + 1), ctx->sbcd, bad_msg[bad_dip]);
      }
      memcpy(&ctx->dipinfo[i], &ctx->dipinfo0, sizeof(SDIP));
      iv2 = iv1 = ctx->glob.vibfac;
      isym = 0;
    }
    isimag[i] = 0;
    ii = iv2 * ifac + iv1;
    ijv = ((unsigned int)ii << 2) + (unsigned int)isym;
    if (ndecv < 3)
    {
      idip[ibcd + 5] = idip[ibcd + ndecv + 2];
      idip[ibcd + 4] = idip[ibcd + ndecv + 1];
      idip[ibcd + 3] = (bcd_t)0;
    }
    else
    {
      idip[ibcd + 4] &= 0xf0;
      idip[ibcd + 3] = (bcd_t)((ijv >> 16) & 0xff);
    }
    idip[ibcd + 2] = (bcd_t)(ijv >> 8);
    idip[ibcd + 1] = (bcd_t)ijv;
    if (si1 == 0)
      idip[ibcd + 4] &= 0x0f;
  }
  kmin = 0;
  nmin = ndip;
  imag = 0;
  for (kd = 0; kd < 8; ++kd)
  {
    ld = kd ^ ctx->glob.stdphase;
    iv1 = 0;
    iv2 = 0;
    nn = 0;
    for (isym = 1; isym <= 3; ++isym)
    {
      i = ctx->ipwr2[isym];
      if (TEST(ctx->glob.phasemask, i) && TEST(kd, i))
        break;
      ii = isym;
      if (TEST(ld, i))
        ii += 4;
      nn += nimag[ii + 8];
      iv1 += nimag[ii ^ 4];
      iv2 += nimag[ii];
    }
    if (isym <= 3)
      continue;
    if (iv1 < iv2)
    {
      nn += iv1;
      ii = 1;
    }
    else
    {
      nn += iv2;
      ii = 0;
    }
    if (nn < nmin)
    {
      nmin = nn;
      kmin = kd;
      imag = ii;
      if (nmin == 0)
        break;
    }
  }
  ctx->glob.stdphase = (kmin ^ ctx->glob.stdphase) & 7;
  if (kmin != 0)
  { /* update phases */
    fprintf(lu,
            "NON-STANDARD PHASE CONVENTION IN USE, %2d\n", ctx->glob.stdphase);
    for (isym = 1; isym <= 3; ++isym)
    {
      ctx->ixphase[isym] = (ctx->ixphase[isym] + kmin) & 1;
      kmin = kmin >> 1;
    }
  }
  if (imag != 0)
  { /* fix flags */
    for (i = 0; i < ndip; ++i)
    {
      iflg = ctx->dipinfo[i].flg;
      if ((iflg & MELEC) == 0)
        continue;
      iflg |= MDIPI;
      ctx->dipinfo[i].flg = (short)iflg;
      isimag[i] = 1;
    }
  }
  if (nmin > 0)
  {
    for (i = 0; i < ndip; ++i)
    {
      iflg = ctx->dipinfo[i].flg;
      if ((iflg & MIMAG) == 0)
        continue;
      isym = (int)(idip[ibcd + 1] & 3);
      if (isym == 0)
        continue;
      ii = ctx->ixphase[isym] + isym;
      if (kmin != 0 && TEST(iflg, MELEC))
        ++ii;
      if (TEST(iflg, MODD))
        ++i;
      if (EVEN(ii))
        continue;
      /* set dipole to zero */
      putbcd(ctx->sbcd, NSBCD, &idip[ibcd]);
      fprintf(lu,
              " WARNING: dipole %3d %s has bad phase and is ignored\n",
              (i + 1), ctx->sbcd);
      ctx->dipinfo[i].fac = 0.;
      ctx->dipinfo[i].flg = (short)0;
      iv1 = ctx->glob.vibfac;
      isym = 0;
      ii = iv1 * ifac + iv1;
      ijv = ((unsigned int)ii << 2) + (unsigned int)isym;
      if (ndecv < 3)
      {
        idip[ibcd + 3] = (bcd_t)0;
      }
      else
      {
        idip[ibcd + 4] &= 0xf0;
        idip[ibcd + 3] = (bcd_t)((ijv >> 16) & 0xff);
      }
      idip[ibcd + 2] = (bcd_t)(ijv >> 8);
      idip[ibcd + 1] = (bcd_t)ijv;
      if (--nmin == 0)
        break;
    }
  }
  return 0;
} /* setint */

/**
 * @brief Reads and sets up global options for calculations from an input file.
 *        This includes vibrational state definitions, spin configurations, K ranges,
 *        statistical weights, symmetry options, and diagonalization preferences.
 * @param input_file_unit File pointer for the input option file.
 * @param num_qn_formats_out Output: Pointer to store the number of quantum number formats determined.
 * @param molecule_type_out Output: Pointer to store molecule type (2 for linear, 3 for non-linear). Renamed itd.
 * @param num_bcd_digits_param_id_out Output: Pointer to store number of BCD digits per parameter ID. Renamed nbcd.
 * @param parameter_name_file_path_out Output: Character array to store the name of the parameter name file (e.g., "sping.nam"). Renamed namfil.
 * @return int Number of option cards read from the input file, or -1 if end-of-file is reached prematurely.
 */
int setopt(struct SpinvContext *ctx, FILE *luin, int *nfmt, int *itd, int *nbcd, char *namfil)
{
  /* this subroutine reads in quantum quantities */
  /*     on entry: */
  /*         LUIN is the fortran unit number to read from */
  /*         NFMT is the maximum number of quanta */
  /*     on return: */
  /*         NFMT= number of quantum number formats for catalog */
  /*         ITD= 2 for linear, 3 for non-linear */
  /*         NAMFIL contains a file name on which to find parameter names */
  /*         SETOPT = number of cards read, -1 for end-of-file */
#define NVEC 11
#define NDECSPIN 11
#define MAXIAX 11
#define NCARD 130
  static int isymv[] = {2, 2, 2, 2, 2, 2, 3, 4, 4, 5, 6, 6};
  static int iaxv[] = {2, 1, 2, 3, 4, 5, 2, 4, 5, 2, 4, 5};
  SVIB *pvinfo, *pvinfom;
  short *pshort;
  double rvec[NVEC], rvec0[NVEC], vsym;
  size_t nl;
  int ivib, lopt, i, ii, j, k, nvibm, isym, iend, ivsym;
  int iwtpl, iwtmn, nt, ncards, knnmin, knnmax, si2, iax;
  int spinsgn, ewt0, ewt1, nspinqn;
  int ivwt[3];
  bcd_t bcdspin[NDECSPIN];
  char card[NCARD], ctyp;

  ctx->cgetv[0].cblk = ctx->cgetv[1].cblk = 0; /* reset store for getqn */
  cjjini();
  ctx->glob.lsym = TRUE;
  ctx->glob.esym = TRUE;
  ctx->glob.esymdec = 100;
  ctx->glob.maxwt = 0;
  ctx->glob.stdphase = 0;
  ctx->glob.newlz = FALSE;
  ctx->glob.nofc = FALSE;
  ctx->glob.g12 = 0;
  rvec0[0] = 1;
  rvec0[1] = 0;
  rvec0[2] = (double)MAXN_DIRCOS;
  rvec0[3] = 0;
  rvec0[4] = 1;
  rvec0[5] = 1;
  rvec0[6] = 1;
  rvec0[7] = 0;
  rvec0[8] = 99;
  rvec0[9] = 0;
  rvec0[10] = 0;
  ncards = nvibm = nt = si2 = nspinqn = 0;
  *itd = 2;
  bcdspin[0] = (bcd_t)NDECSPIN;
  /*     read option cards */
  do
  {
    if (fgetstr(card, NCARD, luin) <= 0)
      return -1;
    ++ncards;
    ctyp = C0;
    if (isalpha((int)card[0]))
      ctyp = card[0];
    else if (isalpha((int)card[1]))
      ctyp = card[1];
    iend = getbcd(card, bcdspin, NDECSPIN);
    spinsgn = NEGBCD(bcdspin[0]);
    dcopy(NVEC, rvec0, 1, rvec, 1);
    rvec[7] = 0.;
    pcard(&card[iend], rvec, NVEC, NULL);
    lopt = (int)rvec[0];
    knnmin = (int)rvec[1];
    knnmax = (int)rvec[2];
    iax = (int)rvec[4];
    iwtpl = (int)rvec[5];
    iwtmn = (int)rvec[6];
    vsym = rvec[7];
    ewt0 = (int)rvec[8];
    ewt1 = 0;
    if (vsym < -0.5)
    {
      ivsym = -1;
      vsym = 0.;
    }
    else if (vsym > 0.5 && vsym < 2.e+15)
    {
      ivsym = 1;
    }
    else
    {
      ivsym = 0;
      vsym = 0.;
    }
    ivib = lopt;
    if (ivib < 0)
      ivib = -ivib;
    if (knnmax < 0 || knnmax > MAXN_DIRCOS)
      knnmax = MAXN_DIRCOS;
    if (knnmin > knnmax)
      knnmax = knnmin;
    if (knnmin != knnmax)
      *itd = 3;
    if (ncards == 1)
    {
      dcopy(NVEC, rvec, 1, rvec0, 1);
      strcpy(namfil, "sping.nam");
      if (ctyp != C0)
        namfil[4] = ctyp;
      ctx->glob.ixz = (int)rvec[3];
      ctx->glob.idiag = (int)rvec[9];
      k = (int)rvec[10];
      ctx->glob.stdphase = MOD(k, 10);
      k = k / 10;
      if (ODD(k))
        ctx->glob.newlz = TRUE;
      if (ODD2(k))
      {
        ctx->glob.nofc = TRUE;
      }
      else if ((k & 4) != 0)
      {
        ctx->glob.g12 = 1;
      }
      if (ivib > 80)
      {
        if (ivib > MAXVIB)
          ivib = MAXVIB;
        if (ivib > 80 && sizeof(int) == 2)
          ivib = 80;
      }
      else if (ivib <= 0)
      {
        ivib = 1;
      }
      if (ctx->glob.nvib > 1)
      {
        free(ctx->vinfo);
        ctx->vinfo = NULL;
      }
      ctx->vinfo1.spt = ctx->sptzero;
      if (ivib > 1)
      {
        nvibm = ivib - 1;
        nl = (size_t)ivib * sizeof(SVIB);
        ctx->vinfo = (SVIB *)mallocq(nl);
        memcpy(ctx->vinfo, &(ctx->vinfo1), sizeof(SVIB));
      }
      else
      {
        ctx->vinfo = &(ctx->vinfo1);
      }
      for (k = 0; k <= 4; ++k)
        ctx->vinfo->wt[k] = 0;
      ctx->vinfo->ewt[0] = ctx->vinfo->ewt[1] = 0;
      ctx->glob.nvib = ivib;
      if (ivib <= 9)
      {
        ctx->glob.vibfac = 9;
        ctx->glob.vibdec = 1;
        k = 4; /* 16 */
      }
      else if (ivib <= 99)
      {
        ctx->glob.vibfac = 99;
        ctx->glob.vibdec = 2;
        k = 7; /* 128 */
      }
      else
      {
        ctx->glob.vibfac = 999;
        ctx->glob.vibdec = 3;
        k = 10; /* 1024 */
      }
      ctx->glob.msshft = (unsigned int)(k + 2);
      /* msshft = shifts needed for vib and symmetry */
      ctx->glob.msmask = (int)(1 << (unsigned int)k) - 1;
      ivib = 0;
      ctx->glob.oblate = (lopt < 0);
      ctx->glob.nqnn = 3;
      if (spinsgn != 0)
        ctx->glob.nqnn = 1;
      if (ewt0 < 0)
        ctx->glob.esymdec = 1000;
    }
    else if (ivib >= ctx->glob.nvib)
    {
      /* if IVIB out of range read another card */
      continue;
    }
    pvinfo = &ctx->vinfo[ivib];
    k = iax;
    if (iax < 0)
      iax = -iax;
    if (iax > MAXIAX)
      iax = 0;
    if (iax <= 3 && ctx->glob.oblate)
      iax = ctx->revsym[iax];
    isym = isymv[iax];
    iax = iaxv[iax];
    pvinfo->gsym = (short)(isym << 1);
    if (isym >= 3)
    {
      j = 2 - (isym & 1);
      if (ctx->glob.maxwt < j)
        ctx->glob.maxwt = j;
    }
    /*  set up minimum K values */
    pvinfo->knmax = (short)knnmax;
    pvinfo->knmin[0] = 0;
    pvinfo->knmin[1] = 1;
    pvinfo->knmin[2] = 1;
    pvinfo->knmin[3] = 2;
    if (k < 0)
    { /* iax was negative */
      if (isym >= 3)
      {
        pvinfo->knmin[0] = -2;
        pvinfo->knmin[3] = -2;
      }
      pvinfo->gsym |= 1;
    }
    if (knnmax != 0 && ctx->glob.nqnn == 1)
      ctx->glob.nqnn = 2;
    /* set ewt */
    if (ewt0 < 0)
      ewt0 = -ewt0;
    if (ewt0 >= ctx->glob.esymdec)
    {
      ewt1 = ewt0 / ctx->glob.esymdec;
      ewt0 -= ctx->glob.esymdec * ewt1;
    }
    if (isym < 3)
    {
      ewt0 = ctx->glob.esymdec - 1;
      ctx->glob.esym = FALSE;
    }
    else if (isym == 4 && iax == 5)
    { /* 4-fold B entry has no E */
      ewt0 = -1;
    }
    if (isym > 3 && pvinfo->wt[4] < 4)
    {
      pvinfo->ewt[0] = pvinfo->ewt[1] = (short)0;
      pvinfo->wt[4] |= 4;
    }
    if (isym == 6)
    { /* six-fold, iax = 4,5 */
      pvinfo->ewt[iax - 4] = (short)ewt0;
    }
    else if (ewt0 >= 0)
    {
      pvinfo->ewt[0] = pvinfo->ewt[1] = (short)ewt0;
    }
    pvinfo->lvupper = 0;
    if (ewt1 != 0)
    {
      if (ewt1 > 10)
      {
        pvinfo->lvupper = 2;
        ewt1 -= 10;
      }
      k = isym;
      if (k < 3)
        k = 3;
      if (ewt1 > k)
        ewt1 = MOD(ewt1 - 1, k) + 1;
      pvinfo->lvupper += (short)(ewt1 << 2);
      pvinfo->knmin[0] = 0;
      pvinfo->knmin[3] = 0; /* 0,1,1,0 */
    }
    pvinfo->lvqn = (short)ewt1;
    for (i = 0; i <= 3; ++i)
    {
      k = pvinfo->knmin[i];
      if (k < 0 && knnmin == 0)
        continue;
      while (k < knnmin)
        k += 2;
      pvinfo->knmin[i] = (short)k;
    }
    k = getsp(ctx, bcdspin, pvinfo);
    if (nspinqn < k)
      nspinqn = k;
    pvinfo->nqn = (short)k;
    if (ncards == 1)
    { /* fill in defaults with wt = 0 */
      pvinfom = ctx->vinfo;
      for (i = 1; i <= nvibm; ++i)
      {
        ++pvinfom;
        memcpy(pvinfom, ctx->vinfo, sizeof(SVIB));
      }
      ivib = -ctx->glob.nvib;
    }
    setwt(pvinfo, ivib, iax, iwtpl, iwtmn, vsym); /* set up weights */
  } while (ivsym < 0); /*   end of reading */
  setsp(ctx);
  /* set up format */
  if (ctx->glob.nqnn == 3)
    *itd = 3;
  ctx->glob.nqn = ctx->glob.nqnn;
  ctx->glob.vibfmt = (nvibm != 0);
  ctx->glob.iqfmt0 = 0;
  if (ctx->glob.vibfmt)
  {
    ctx->glob.nqn += 1;
    ctx->glob.iqfmt0 = 1000;
  }
  nt = ctx->glob.nqn;
  k = nt + nspinqn;
  ctx->glob.nqn = k;
  if (k > (*nfmt))
    k = nt + 2;
  ctx->glob.maxqn = k;
  ctx->glob.nqn0 = nt;
  ctx->glob.iqfmt0 += 100 * nt + k;
  /*     make sure +/- l values come in pairs */
  pvinfo = pvinfom = ctx->vinfo;
  k = 0;
  for (i = 0; i <= nvibm; ++i)
  {
    if (i > 0)
    {
      pvinfom = pvinfo;
      ++pvinfo;
    }
    pvinfo->nqn += (short)nt;
    pvinfo->wt[4] &= 3;
    j = pvinfo->wt[4];
    if (EVEN(pvinfo->gsym))
    {
      pvinfo->wt[4] = 0; /* no Itot */
    }
    else if (pvinfo->wt[0] == pvinfo->wt[j] &&
             pvinfo->wt[1] == pvinfo->wt[j ^ 1])
    {
      pvinfo->wt[4] = (short)((pvinfo->gsym > 6) ? -1 : 0);
    }
    j = pvinfo->lvqn;
    if (k != 0)
    { /* previous l was non-zero and first of pair */
      isym = (int)((unsigned int)pvinfo->gsym >> 1);
      if (isym >= 3)
      {
        if ((k + k) == isym && j != k)
        {        /* B symmetry */
          k = j; /* single instance */
          pvinfom->lvqn = 0;
          pshort = pvinfom->wt;
          ii = pshort[0];
          pshort[0] = pshort[2];
          pshort[2] = (short)ii;
          ii = pshort[1];
          pshort[1] = pshort[3];
          pshort[3] = (short)ii;
          pshort = pvinfom->ewt;
          ii = pshort[0];
          pshort[0] = pshort[1];
          pshort[1] = (short)ii;
          continue;
        }
        else
        {
          if (MOD(j + k, isym) != 0)
            break;
        }
      }
      if (k > j)
      {
        if (j == 0)
          break;
        pvinfom->lvqn = (short)(-j);
        if (pvinfo->knmin[0] == 0)
          pvinfo->knmin[0] = 2;
        if (pvinfo->knmin[3] == 0)
          pvinfo->knmin[3] = 2;
      }
      else
      {
        pvinfo->lvqn = (short)(-k);
        if (pvinfom->knmin[0] == 0)
          pvinfom->knmin[0] = 2;
        if (pvinfom->knmin[3] == 0)
          pvinfom->knmin[3] = 2;
      }
      for (j = 0; j <= 4; ++j)
      {
        if (pvinfo->wt[j] != pvinfom->wt[j])
          break;
      }
      if (j <= 4)
        break;
      if (pvinfo->gsym != pvinfom->gsym)
        break;
      if (pvinfo->ewt[0] != pvinfom->ewt[0])
        break;
      if (pvinfo->ewt[1] != pvinfom->ewt[1])
        break;
      pvinfo->lvupper += 1;
      j = 0;
    }
    k = j;
  }
  if (k != 0)
  {
    puts("L DOUBLETS SHOULD BE IN PAIRS");
    exit(EXIT_FAILURE);
  }
  k = 0;
  pvinfo = ctx->vinfo;
  for (ivib = 0; ivib <= nvibm; ++ivib)
  {
    nt = pvinfo->nspstat;
    ii = setgsym(ctx, (int)pvinfo->gsym);
    for (isym = 0; isym < 4; ++isym)
    {
      if (getwt(ctx, pvinfo, isym, 0, ivwt) <= 0)
        continue;
      for (ii = 1; ii <= nt; ++ii)
      {
        if (getwt(ctx, pvinfo, isym, ii, ivwt) > 0)
          ++k;
      }
    }
    ++pvinfo;
  }
  if (ctx->glob.nbkpj > 0)
  {
    free(ctx->moldv);
    ctx->moldv = NULL;
    free(ctx->blkptr);
    ctx->blkptr = NULL;
  }
  ctx->glob.nbkpj = k;
  nl = (size_t)k * sizeof(int);
  ctx->moldv = (int *)mallocq(nl);
  ctx->moldv[0] = 0;
  nl += sizeof(int);
  ctx->blkptr = (int *)mallocq(nl);
  ctx->blkptr[0] = 0;
  k = 0;
  pvinfo = ctx->vinfo;
  for (ivib = 0; ivib <= nvibm; ++ivib)
  {
    nt = pvinfo->nspstat;
    ii = setgsym(ctx, (int)pvinfo->gsym);
    for (isym = 0; isym < 4; ++isym)
    {
      if (getwt(ctx, pvinfo, isym, 0, ivwt) <= 0)
        continue;
      for (ii = 1; ii <= nt; ++ii)
      {
        if (getwt(ctx, pvinfo, isym, ii, ivwt) > 0)
        {
          ctx->blkptr[k] = k;
          ctx->moldv[k] = isym + (ivib << 2) + (ii << ctx->glob.msshft);
          ++k;
        }
      }
    }
    ++pvinfo;
  }
  ctx->blkptr[k] = k;
  pvinfo = NULL;
  *nbcd = (NDECPAR + 1) + ctx->glob.vibdec;
  if (sizeof(int) == 2 && ctx->glob.vibdec > 2)
  {
    puts("too many vibrations for 16-bit integer computer");
    ncards = 0;
  }
  *nfmt = 1;
  if (nspinqn != 0 && nvibm != 0)
    *nfmt = setfmt(ctx, &k, -1);
  return ncards;
} /* setopt */

/**
 * @brief Determines the quantum number output format code (QNFMT) for vibrational states.
 *        If nfmt_or_vib_limit > 0, it fills iqnfmt_out array with format codes for each vibrational state.
 *        If nfmt_or_vib_limit = -1, it checks if all states up to glob.nvib have the same format code
 *        and returns 1 if they do, or glob.nvib if they differ (this signals multiple formats needed).
 * @param iqnfmt_out Output: If nfmt_or_vib_limit > 0, this array is filled with QNFMT codes for each vibrational state.
 *                   If nfmt_or_vib_limit = -1, this parameter is effectively a dummy (value passed as &k from setopt).
 * @param nfmt_or_vib_limit If > 0, it's the limit of vibrational states to process for filling iqnfmt_out.
 *                          If -1, it's a flag to check format consistency across all vibrational states.
 * @return int If nfmt_or_vib_limit > 0, returns glob.maxqn (total QN fields to be output).
 *             If nfmt_or_vib_limit = -1, returns 1 if all formats are same, or a value > 1 (e.g. num_vib_states_to_check) if they differ.
 */
int setfmt(struct SpinvContext *ctx, int *iqnfmt, int nfmt)
{
  SVIB *pvinfo;
  short *jjs;
  int iqfmtq, iqfmt0, i, j, k, iv, si2, ff, nset, nspinv, nvib;
  iqfmt0 = 0;
  nvib = ctx->glob.nvib;
  if (nfmt > 0 && nfmt < nvib)
    nvib = nfmt;
  pvinfo = ctx->vinfo;
  for (iv = 0; iv < nvib; ++iv)
  {
    iqfmtq = ctx->glob.iqfmt0;
    i = setgsym(ctx, (int)pvinfo->gsym);
    jjs = pvinfo->spt;
    nset = jjs[0];
    nspinv = nset - 1;
    if (nspinv > ctx->nspin)
      nspinv = ctx->nspin;
    jjs += nset;
    ff = 200 - jjs[0];
    if ((int)pvinfo->nqn > ctx->glob.maxqn)
    {
      iqfmtq += 10 * (ff & 1) + 4000;
    }
    else
    {
      si2 = 0;
      k = ctx->glob.maxqn - ctx->glob.nqn0 - 1;
      for (i = 0; i < k; ++i)
      {
        if (i == ctx->itsym && ctx->itsym < ctx->itptr)
        {
          i = ctx->itptr;
          j = 0;
          si2 = si2 << 1;
        }
        if (i < nspinv)
        {
          j = jjs[i + 1];
          if (i < ctx->itsym)
            j += ff;
        }
        else
        {
          j = ff;
        }
        si2 = (si2 + (j & 1)) << 1;
      }
      si2 += (ff & 1);
      if (si2 > 9)
        si2 &= 7;
      if (ctx->itptr < ctx->nspin)
        iqfmtq += 2000;
      if (ctx->itsym < ctx->itptr)
        iqfmtq += 4000;
      iqfmtq += 10 * si2;
    }
    if (iv == 0)
      iqfmt0 = iqfmtq;
    if (nfmt < 0)
    {
      if (iqfmtq != iqfmt0)
        return nvib;
    }
    else
    {
      iqnfmt[iv] = iqfmtq;
    }
    ++pvinfo;
  }
  if (nfmt <= 0)
    return 1;
  return ctx->glob.maxqn;
} /* setfmt */

/**
 * @brief Initializes global parameters and allocates/frees memory related to Hamiltonian blocks
 *        and parameter processing. Called from main program to set up before calculations or
 *        to clean up.
 * @param list_unit File pointer for listing output. Renamed lu.
 * @param num_total_params Total number of parameters. Renamed npar.
 * @param bcd_parameter_ids Array of BCD parameter identifiers. Renamed idpar.
 * @param parameter_values Array of parameter values. Renamed par.
 * @param max_f_qn_in_block Pointer to store maximum F quantum number used in any block. Renamed nbkpf.
 * @param max_hamiltonian_dim Pointer to store maximum Hamiltonian dimension. Renamed negy.
 * @return int Returns a value related to vibrational state packing (vibfac^2 related), or 0 if npar is 0.
 *         The exact meaning of the return value might be specific to the calling context.
 */
int setblk(struct SpinvContext *ctx, FILE *lu, const int npar, bcd_t *idpar, const double *par, int *nbkpf, int *negy)
{
  /* this subroutine initializes  quantum quantities */
  /*     on entry: */
  /*         LU= logical unit for listing */
  /*         IDPAR= parameter id */
  /*         NBKPF= maximum value of F */
  /*         NEGY = max dimensioned energy in main program */
  /*     on return: */
  /*         NBKPF= blocks per F */
  /*         NEGY= possibly downsized NEGY to account for local storage */
  static char csym[] = {'x', 'c', 'b', 'a'};
  SVIB *pvinfo;
  SPAR *spar_now;
  short *jjs;
  double ai;
  size_t nl;
  int nodd, isym, jsym, ixzt, i, k, kl, lt, n, kdiag, nvibm, neven, nsize, alpha;
  int kd, ld, ii, jj, im, kk, ni, nj, iv, jv, ivmin, ix, lv, lx, npair, ldel;
  int knnmin, ntstat, si1, si2, iff, ikq, jjt, ins, neuler, sznz, iff0;
  int iwt, gsym, maxblk, maxsblk, itmp, nd, ifc, nf, kavg;
  unsigned int ivsym;
  int ivwt[3], jvwt[3];
  char ctmp;

  ctx->cgetv[0].cblk = ctx->cgetv[1].cblk = 0; /* reset store for getqn */
  kdiag = ctx->glob.idiag;
  if (npar <= 0)
  {
    if (ctx->glob.nvib > 1)
      free(ctx->vinfo);
    ctx->vinfo = &ctx->vinfo1;
    ctx->glob.nvib = 0;
    if (ctx->glob.nbkpj > 0)
    {
      free(ctx->moldv);
      ctx->moldv = &ctx->zmoldv;
      free(ctx->blkptr);
      ctx->blkptr = &ctx->zblkptr;
      ctx->glob.nbkpj = 0;
    }
    if (ctx->glob.maxblk > 0)
    {
      free(ctx->ivs);
      ctx->ivs = &ctx->zivs;
      free(ctx->ikmin);
      ctx->ikmin = &ctx->zikmin;
      free(ctx->ibkptr);
      ctx->ibkptr = &ctx->zibkptr;
      ctx->glob.maxblk = 0;
    }
    if (ctx->ndmx > 0)
    {
      free(ctx->iqnsep);
      ctx->iqnsep = &ctx->ziqnsep;
      free(ctx->idx);
      ctx->idx = &ctx->zidx;
      free(ctx->jdx);
      ctx->jdx = &ctx->zjdx;
      free(ctx->wk);
      ctx->wk = &ctx->zwk;
      ctx->ndmx = 0;
    }
    return 0;
  }
  nvibm = ctx->glob.nvib - 1;
  n = (*nbkpf);
  if (n > 0)
  {
    iff0 = n << 1;
    pvinfo = ctx->vinfo;
    for (i = 0; i <= nvibm; ++i)
    {
      jjs = pvinfo->spt;
      k = jjs[0];
      /*  select largest N relative to F */
      ni = (iff0 - jjs[k]) >> 1;
      if ((int)pvinfo->knmax > ni)
        pvinfo->knmax = (short)ni;
      ++pvinfo;
    }
    pvinfo = NULL;
  }
  if (pasort(ctx, lu, npar, idpar, par) == 0)
    ctx->glob.idiag = -1;
  pvinfo = NULL;
  iff0 = 200;
  /*    find interactions */
  ntstat = ctx->glob.nbkpj;
  sznz = ctx->ndmax = 0;
  if (ntstat <= 0)
    return 0;
  nf = ctx->nspin - 1;
  for (ii = 0; ii < ntstat; ++ii)
  {
    im = ctx->moldv[ii];
    iff = iff0;
    kk = getqs(ctx, im, iff, -1, -1, ctx->ixcom, ctx->iscom, &iv);
    if (ODD(ctx->iscom[ctx->nspin]))
    {
      /* special for odd spin multiplicity (make N integer) */
      kk = getqs(ctx, im, --iff, 0, -1, ctx->ixcom, ctx->iscom, &iv);
    }
    isym = ctx->ixcom[XSYM];
    ni = ctx->ixcom[XNVAL];
    pvinfo = &ctx->vinfo[ctx->ixcom[XVIB]];
    gsym = pvinfo->gsym;
    /* nsym = */ setgsym(ctx, gsym);
    getwt(ctx, pvinfo, isym, kk, ivwt);
    for (jj = 0; jj < ii; ++jj)
    {
      if (ctx->blkptr[ii] == ctx->blkptr[jj])
        continue;
      jjt = ctx->moldv[jj];
      kk = getqs(ctx, jjt, iff, 0, -1, ctx->jxcom, ctx->jscom, &jv);
      jsym = ctx->jxcom[XSYM];
      nj = ctx->jxcom[XNVAL];
      nd = ni - nj;
      pvinfo = &ctx->vinfo[ctx->jxcom[XVIB]];
      if (pvinfo->gsym != (short)gsym)
        continue;
      getwt(ctx, pvinfo, jsym, kk, jvwt);
      for (k = ctx->glob.maxwt; k >= 0; --k)
      {
        if (ivwt[k] != jvwt[k])
          break;
      }
      if (k >= 0)
        continue;
      /* check for connection restrictions on N */
      ixzt = ctx->glob.ixz;
      if (ODD(ixzt) && ctx->iscom[ctx->nspin] != ctx->jscom[ctx->nspin])
        continue;
      for (k = 0; k < nf; ++k)
      {
        if (k == ctx->itsym)
          k = ctx->itptr;
        ixzt = ixzt >> 1;
        /* check for matching spin multiplicity */
        if (ODD(ctx->iscom[k] + ctx->jscom[k]))
          break;
        /* check for connection restrictions on spin */
        if (ODD(ixzt) && ctx->iscom[k] != ctx->jscom[k])
          break;
      }
      if (k < nf)
        continue;
      ivmin = (iv < jv) ? iv : jv;
      itmp = iv + jv + ivmin * ctx->glob.vibfac;
      ivsym = blksym(ctx->ixcom, ctx->jxcom) + ((unsigned int)itmp << 2);
      for (spar_now = ctx->spar_head[ivmin]; spar_now != NULL;
           spar_now = spar_now->next)
      {
        if (spar_now->ipsym < ivsym)
          break;
        if (spar_now->ipsym > ivsym)
          continue;
        if (TEST(spar_now->flags, MNSQ))
          continue;
        if (sznz < 0)
          sznzfix(ctx, sznz, ni, nj, ctx->ixcom, ctx->jxcom, ctx->iscom, ctx->jscom);
        kl = idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins, &si1, &si2,
                    &sznz, &ifc, &alpha, &ldel, &kavg);
        if (ODD(neuler))
          continue;
        if (sznz > 0)
        {
          sznz = sznzfix(ctx, sznz, ni, nj, ctx->ixcom, ctx->jxcom, ctx->iscom, ctx->jscom);
        }
        if (getmask(ctx, ctx->ixcom, ctx->jxcom, kd, ldel, kl, alpha) == 0)
          continue;
        npair = getll(ctx, 0, ld, lt, kd, si1, si2, ctx->lscom, ctx->iscom, ctx->jscom);
        if (npair < 0)
          continue;
        if (kd > 0)
          ctx->glob.idiag = kdiag;
        if (nd > ctx->ndmax)
          ctx->ndmax = nd;
        kd = ctx->blkptr[ii] - ctx->blkptr[jj];
        if (kd != 0)
        { /* II and JJ are coupled */
          if (kd > 0)
          { /* rename ptr II to ptr JJ */
            kk = ctx->blkptr[ii];
            ix = jj;
          }
          else
          { /* rename ptr JJ to ptr II */
            kk = ctx->blkptr[jj];
            ix = ii;
          }
          ctx->glob.idiag = kdiag;
          for (k = 0; k <= ii; ++k)
          {
            if (ctx->blkptr[k] == kk)
              ctx->blkptr[k] = ctx->blkptr[ix];
          }
        }
        break;
      } /* loop over parameters */
    } /* jj loop */
  } /* ii loop */
  /*  block connections now found */
  if (ctx->glob.idiag < 0)
    ctx->glob.idiag = -1;
  switch (ctx->glob.idiag)
  {
  case -1:
    fputs("NO DIAGONALIZATION NEEDED\n", lu);
    break;
  default:
    fputs("ENERGY SORT OF WANG SUB-BLOCKS\n", lu);
    break;
  case 1:
    fputs("EIGENVECTOR SORT OF STATES\n", lu);
    break;
  case 2:
    fputs("ENERGY FOLLOWS FIRST ORDER WITHIN WANG SUB-BLOCK\n", lu);
    break;
  case 3:
    fputs("ENERGY FOLLOWS ORDER OF TAU WITHIN V AND SPIN SUB-BLOCKS\n", lu);
    break;
  case 4:
    fputs("ENERGY FOLLOWS K*K WITHIN V AND SPIN SUB-BLOCKS\n", lu);
    break;
  case 5:
    fputs("ENERGY FOLLOWS FIRST ORDER WITHIN V AND SPIN SUB-BLOCKS\n", lu);
  }
  if (ctx->glob.stdphase > 0)
    fprintf(lu, "NON-STANDARD PHASE CONVENTION IN USE, %2d\n", ctx->glob.stdphase);
  if (ctx->glob.newlz)
    fputs("Lz DEFINED EXPLICITLY\n", lu);
  if (ctx->glob.nofc)
    fputs("ALTERNATE DEFINITION OF PARAMETER FC FIELD\n", lu);
  if (ctx->glob.g12 != 0)
    fputs("G12 group alternating sign of Fourier coefficients with K  \n", lu);
  if (ctx->glob.oblate)
  {
    fputs("OBLATE ROTOR\n", lu);
    ctmp = csym[1];
    csym[1] = csym[3];
    csym[3] = ctmp;
  }
  else
  {
    fputs("PROLATE ROTOR\n", lu);
  }
  if (ctx->glob.nqnn == 1)
  {
    fputs("LINEAR MOLECULE QUANTA, K SUPPRESSED\n", lu);
  }
  else if (ctx->glob.nqnn == 2)
  {
    fputs("SYMMETRIC TOP QUANTA\n", lu);
  }
  fputs("    V KMIN KMAX WTPL WTMN ESYMWT NSYM SPINS\n", lu);
  pvinfo = ctx->vinfo;
  jj = 0;
  for (iv = 0; iv <= nvibm; ++iv)
  {
    knnmin = pvinfo->knmin[0];
    ii = pvinfo->knmin[1];
    if (ii < knnmin)
      knnmin = ii;
    ii = pvinfo->knmin[2];
    if (ii < knnmin)
      knnmin = ii;
    ii = pvinfo->knmin[3];
    if (ii < knnmin)
      knnmin = ii;
    if (knnmin < 0)
      knnmin = 0;
    k = pvinfo->gsym;
    isym = k >> 1;
    ii = (int)(pvinfo->lvupper >> 2);
    if ((pvinfo->lvupper & 2) != 0)
      ii += 10;
    if (ODD(k))
    {
      isym = -isym;
      jj = 1;
    }
    ii = ii * ctx->glob.esymdec + pvinfo->ewt[0];
    fprintf(lu, " %4d %4d %4d %4d %4d %6d %4d", iv, knnmin, pvinfo->knmax,
            pvinfo->wt[0], pvinfo->wt[3], ii, isym);
    jjs = pvinfo->spt;
    ni = jjs[0];
    for (i = 1; i < ni; ++i)
    {
      ai = jjs[i] * 0.5;
      fprintf(lu, "%5.1f", ai);
    }
    fputc('\n', lu);
    if (EVEN(isym))
    { /* even */
      i = ii;
      ii -= pvinfo->ewt[0];
      if (isym != 4)
        ii += pvinfo->ewt[1];
      if (i != ii || pvinfo->wt[0] != pvinfo->wt[2] ||
          pvinfo->wt[3] != pvinfo->wt[1])
      {
        fprintf(lu, " %4d %4d %4d %4d %4d %6d    B\n", iv, knnmin, pvinfo->knmax,
                pvinfo->wt[2], pvinfo->wt[1], ii);
      }
    }
    ++pvinfo;
  }
  if (jj != 0)
  {
    fprintf(lu, "I(TOT) IS LAST QUANTUM BEFORE F AND %s\n",
            "IS INDICATED BY NEGATIVE VALUE OF NSYM ABOVE");
  }
  /* bubble sort block number */
  lx = ntstat;
  while (lx > 1)
  {
    lv = 0;
    im = 0;
    for (i = 1; i < lx; ++i)
    {
      if (ctx->blkptr[im] > ctx->blkptr[i] ||
          (ctx->blkptr[im] == ctx->blkptr[i] && ctx->moldv[im] > ctx->moldv[i]))
      {
        lv = i;
        jj = ctx->blkptr[i];
        ctx->blkptr[i] = ctx->blkptr[im];
        ctx->blkptr[im] = jj;
        jj = ctx->moldv[i];
        ctx->moldv[i] = ctx->moldv[im];
        ctx->moldv[im] = jj;
      }
      im = i;
    }
    lx = lv;
  }
  fputs("BLOCK - WT - SYM - V - TSP - N  -  other quanta  (rel. to F=0 )\n",
        lu);
  /* convert block label to pointer */
  i = ctx->blkptr[0];
  ctx->blkptr[0] = 0;
  n = nsize = neven = nodd = maxblk = maxsblk = 0;
  for (jj = 0; jj <= ntstat; ++jj)
  {
    if (ctx->blkptr[jj] != i)
    {
      ii = jj - ctx->blkptr[n++];
      if (maxblk < ii)
        maxblk = ii;
      ctx->blkptr[n] = jj;
      if (neven > nsize)
        nsize = neven;
      if (nodd > nsize)
        nsize = nodd;
      if (jj == ntstat)
        break;
      i = ctx->blkptr[jj];
      neven = 0;
      nodd = 0;
    }
    jjt = ctx->moldv[jj];
    jjt = getqs(ctx, jjt, 0, -1, 0, ctx->ixcom, ctx->iscom, &ii);
    isym = ctx->ixcom[XSYM];
    iv = ctx->ixcom[XVIB];
    pvinfo = &ctx->vinfo[iv];
    iwt = getwt(ctx, pvinfo, isym, jjt, ivwt);
    kk = pvinfo->knmin[isym];
    if (kk < 0)
    {
      ii = ctx->ixcom[XISYM]; /* get spin symmetry */
      kk = isym;
      if (ctx->is_esym[ii] < 0)
        kk += 2;
      kk &= 2;
    }
    ii = pvinfo->knmax - kk + 2;
    if (ii < 0)
      ii = 0;
    ii >>= 1; /* size of even block */
    neven += ii;
    if (ii > maxsblk)
      maxsblk = ii;
    kk = pvinfo->knmin[3 - isym];
    if (kk < 0)
    {
      ii = ctx->ixcom[XISYM]; /* get spin symmetry */
      kk = isym;
      if (ctx->is_esym[ii] >= 0)
        kk += 2;
      kk &= 2;
    }
    ii = pvinfo->knmax - kk + 2;
    if (ii < 0)
      ii = 0;
    ii >>= 1; /* size of odd block */
    nodd += ii;
    if (ii > maxsblk)
      maxsblk = ii;
    --jjt;
    fprintf(lu, "% 5d %4d    %c %4d %4d", n + 1, iwt, csym[isym], iv, jjt);
    ni = ctx->nspin;
    for (ii = 0; ii < ctx->nspin; ++ii)
    {
      ai = ctx->iscom[ni] * 0.5;
      fprintf(lu, " %5.1f", ai);
      if (ni == ctx->itsym && ctx->itsym < ctx->itptr)
        ii = ctx->itptr;
      ni = ii;
    }
    fputc('\n', lu);
  }
  ++maxblk;
  ctx->glob.maxblk = maxblk;
  fprintf(lu, "Maximum Dimension for Hamiltonian = %d\n", nsize);
  fflush(lu);
  *nbkpf = n;
  ctx->glob.nbkpj = n;
  if (nsize < (*negy))
    *negy = nsize;
  if (ctx->ndmx > 0)
  {
    free(ctx->ivs);
    ctx->ivs = NULL;
    free(ctx->ikmin);
    ctx->ikmin = NULL;
    free(ctx->ibkptr);
    ctx->ibkptr = NULL;
    free(ctx->iqnsep);
    ctx->iqnsep = NULL;
    free(ctx->idx);
    ctx->idx = NULL;
    free(ctx->jdx);
    ctx->jdx = NULL;
    free(ctx->wk);
    ctx->wk = NULL;
    ctx->ndmx = 0;
  }
  ii = maxsblk + maxsblk - 1;
  ctx->ndmx = nsize;
  if (ctx->ndmx < ii)
    ctx->ndmx = ii;
  nl = (size_t)(nsize + ctx->ndmx) * sizeof(double);
  ctx->wk = (double *)mallocq(nl);
  ctx->wk[0] = 0.0;
  /* allocate space for sub-block information */
  nl = (size_t)(ctx->ndmx + 1) * sizeof(short);
  ctx->idx = (short *)mallocq(nl);
  ctx->idx[0] = 0;
  ctx->jdx = (short *)mallocq(nl);
  ctx->jdx[0] = 0;
  nl = nsize * sizeof(short);
  ctx->iqnsep = (short *)mallocq(nl);
  ctx->iqnsep[0] = 0;
  maxblk = maxblk + maxblk;
  nl = maxblk * sizeof(short);
  ctx->ibkptr = (short *)mallocq(nl);
  ctx->ibkptr[0] = 0;
  ctx->ikmin = (short *)mallocq(nl);
  ctx->ikmin[0] = 0;
  nl = (size_t)maxblk * sizeof(int);
  ctx->ivs = (int *)mallocq(nl);
  ctx->ivs[0] = 0;
  ii = ctx->glob.vibfac + 1;
  return (ii * ii);
} /* setblk */

/**
 * @brief Initializes spin factors and sorts/pre-processes Hamiltonian parameters.
 *        This function parses BCD parameter identifiers, determines operator types,
 *        calculates spin-independent and K-independent factors (zfac), and stores
 *        processed parameter information in linked lists (spar_head) for efficient
 *        Hamiltonian construction. It also handles phase conventions for operators.
 * @param list_unit File pointer for listing warnings and information. Renamed lu.
 * @param num_total_params Total number of parameters. Renamed npar.
 * @param bcd_parameter_ids Array of BCD parameter identifiers. Renamed idpar.
 * @param parameter_values Array of parameter values. Renamed par.
 * @return int Returns 1 if any diagonal K-changing operator was found (iret=1), 0 otherwise.
 *         A return of 0 might imply that glob.idiag should be set to -1 (no diagonalization).
 */
int pasort(struct SpinvContext *ctx, FILE *lu, const int npar, bcd_t *idpar, const double *par)
{
  /*     initialize spin factors and arrange operators */
  /*     on entry: */
  /*         NSIZE= block size and dimension */
  /*         IDPAR= list of parameter identifiers ( element 0 is length ) */
  static char *strej[] = {"\n", "(1) vib symmetry does not match\n",
                          "(2) bad spins\n", "(3) bad spin 1 for SzNz\n",
                          "(4) specified K is beyond Kmax\n", "(5) bad weights for Itot\n",
                          "(6) specify prameters only for l >= 0\n",
                          "(7) matrix element MOD(Delta K - Delta l, nsym) is not 0\n",
                          "(8) problem with Itot symmetry\n",
                          "(9) matrix element couples states with differernt weights\n"};

  static int idmy[NDXCOM], initl;
  static bcd_t idunit[NDECPAR];
  SVIB *pvib1, *pvib2;
  SPAR *spar_now, *spar_last, *spar_match, *spar_free;
  short *iiv1, *iiv2;
  double dtmp;
  size_t nl;
  int ifac, iv1d, ipar, iv12q, ikqp, insp, njqt, ityi, isym, neulerp, i, lt;
  int kk, ltp, si1, si2, iv1, iv2, lv1, lv2, kdp, ldp, iv12, ikq, alpha, gsym;
  int ins, neuler, sznz, sznzp, si1p, si2p, nt, kd, ld, ivdif, ibcd, ldel;
  int iret, idflags, ifc, ifcp, ilim, kl, klp, nitot, ndecv, nbcd, imag;
  int alphap, ldelp, kavg, kavgp, notused, nimag[8];
  unsigned int ivsym;
  bcd_t *idval;
  BOOL first;
  const double zero = 0.0;
  const short szero = 0;

  if (ctx->glob.parinit > 0)
  {
    for (i = 0; i < ctx->glob.parinit; ++i)
    {
      spar_free = ctx->spar_head[i];
      while (spar_free != NULL)
      {
        spar_now = spar_free->next;
        free(spar_free);
        spar_free = spar_now;
      }
      ctx->spar_head[i] = spar_free;
    }
    spar_free = spar_now = NULL;
    ctx->glob.parinit = 0;
  }
  if (npar <= 0)
  {
    free(ctx->ipder);
    ctx->ipder = &ctx->zipder;
    initl = 0;
    return 0;
  }
  /*  find reduced spin matrix elements */
  ctx->spfac[0] = 0.;
  ctx->spfac2[0] = 0.;
  ctx->spfac[1] = sqrt(1.5);
  ctx->spfac2[0] = 0.;
  nt = ctx->glob.mxspin;
  for (i = 2; i <= nt; ++i)
  {
    dtmp = 0.5 * i; /* dtmp = I */
    ctx->spfac[i] = sqrt(dtmp * (dtmp + 1.) * (i + 1));
    ctx->spfac2[i] = 0.25 * sqrt((dtmp + 1.5) / (dtmp - 0.5)) / dtmp;
  }
  /* initialize SPECFC and DIRCOS */
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wdiscarded-qualifiers"
  // zero is used as non-const in specfc & dircos but should not actually be modified with these sets of parameters
  specfc(ctx, 0, 0, 0, 0, 0, 0, 0, &zero, &szero, &szero);
  kk = 0;
  dircos(ctx, idmy, idmy, 0, 0, 0, &zero, &szero, &szero, 0, MODD, 0, &kk);
  #pragma GCC diagnostic pop
  ifac = ctx->glob.vibfac + 1;
  ndecv = ctx->glob.vibdec;
  ilim = ifac * ifac - 1;
  if (initl > 0)
  {
    free(ctx->ipder);
    ctx->ipder = NULL;
  }
  nl = (size_t)npar * sizeof(int);
  ctx->ipder = (int *)mallocq(nl);
  initl = npar;
  kk = 1;
  ityi = 1;
  ctx->ipder[0] = 0;
  ibcd = 0;
  nbcd = (int)idpar[0] & 0x7f;
  for (i = 1; i < npar; ++i)
  {
    ibcd += nbcd;
    if (NEGBCD(idpar[ibcd]) == 0)
    {
      ctx->ipder[i] = kk++;
      ityi = i + 1;
    }
    else
    {
      ctx->ipder[i] = -ityi;
    }
  }
  for (i = 0; i < 8; ++i)
    nimag[i] = 0;
  ctx->glob.nfit = kk;
  /* .. code parameter id  for power of N*(N+1) and equivalence */
  /* .. NJQ = POWER + 1  (first occurrance of power) */
  /* .. IPTYP= ITY + 1    (first occurrance of cosine) */
  /* ..      ITY = sequence number for hybrid operators */
  /* .. negate NJQ and IPTYP for successive occurrance */
  /* .. IP points to parameter */
  ctx->glob.parinit = ctx->glob.nvib;
  for (i = 0; i < ctx->glob.nvib; ++i)
  {
    ctx->spar_head[i] = NULL;
  }
  ctx->nsqmax = ityi = iv1d = iret = iv1 = iv2 = nitot = gsym = 0;
  ipar = npar - 1;
  ibcd = ipar * nbcd;
  spar_free = NULL;
  pvib1 = pvib2 = ctx->vinfo;
  while (ipar >= -1)
  {
    notused = -1;
    ityi = 0;
    first = TRUE;
    do
    { /* repeat until parameter subtypes exhausted */
      ivdif = 0;
      if (ipar < 0)
      { /* make sure there is a unit operator */
        idval = idunit;
        iv12q = ilim;
      }
      else
      {
        idval = &idpar[ibcd + 1];
        iv12q = bcd2i(idval[0]);
        if (ndecv > 1)
          iv12q += bcd2i(idval[1]) * 100;
        if (ndecv > 2)
          iv12q += bcd2i(idval[2]) * 10000;
        idval += ndecv;
        iv2 = iv12q / ifac;
        iv1 = iv12q - iv2 * ifac;
        ivdif = iv1 - iv2;
        if (ivdif < 0)
        {
          i = iv1;
          iv1 = iv2;
          iv2 = i;
          iv12q = iv2 * ifac + iv1;
        }
      }
      iv12 = iv12q;
      if (iv12q == ilim)
      {
        iv1 = iv2 = iv1d;
        ivdif = 0;
        iv12 = iv1d * (ifac + 1);
      }
      else if (iv1 >= ctx->glob.nvib)
      {
        /* dummy with v1 = 99 and v2 < 99 */
        if (iv1 == ctx->glob.vibfac)
          notused = 0;
        break;
      }
      if (idval[1] == (bcd_t)0x91 && idval[0] == (bcd_t)0)
      {
        /* special for rho */
        if (iv1 == iv2)
        {
          notused = 0;
          dtmp = par[ipar];
          specfc(ctx, 0, iv1, iv1, 0, 0, 0, 1, &dtmp, &szero, &szero);
        }
        break;
      }
      else if (idval[1] == (bcd_t)0x99 && idval[0] == (bcd_t)0)
      {
        /* dummy parameter = 9900vv' */
        notused = 0;
        break;
      }
      pvib1 = &ctx->vinfo[iv1];
      pvib2 = &ctx->vinfo[iv2];
      gsym = pvib1->gsym;
      if ((short)gsym != pvib2->gsym)
      {
        notused = 1;
        break;
      }
      gsym = setgsym(ctx, gsym);
      nitot = ctx->glob.nitot;
      if (spar_free == NULL)
      { /* allocate more structures */
        spar_free = (SPAR *)mallocq(sizeof(SPAR));
      }
      spar_now = spar_free;
      spar_now->flags = 0;
      spar_now->alpha = C0;
      spar_now->mldel = C0;
      ityi = idpari(ctx, idval, ityi, spar_now);
      if (ityi == 0)
        break;
      njqt = (int)spar_now->njq;
      if (njqt > ctx->nsqmax)
        ctx->nsqmax = njqt;
      isym = (int)spar_now->ipsym;
      idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins, &si1, &si2,
             &sznz, &ifc, &alpha, &ldel, &kavg);
      iiv1 = pvib1->spt;
      iiv2 = pvib2->spt;
      if (si1 > 0 &&
          checksp(ctx, first, si1, si2, iiv1, iiv2, &spar_now->zfac) != 0)
      {
        notused = 2;
        break;
      }
      if (sznz > 0)
      { /*  setup for SzNz */
        if (checksp(ctx, first, 1, 0, iiv1, iiv2, &spar_now->zfac) != 0)
        {
          notused = 3;
          break;
        }
        spar_now->zfac *= 0.5;
      }
      first = FALSE;
      /* factor of 2 for B-C, etc */
      if (kd != 0 && isym == 0)
        spar_now->zfac *= 2.;
      /* setup flags */
      idflags = 1;
      imag = (kd > ld) ? kd : ld;
      if (ifc < 0)
      {
        ++imag;
        isym = 3 - isym;
        idflags |= MFCM;
      }
      if (kavg > pvib1->knmax && kavg > pvib2->knmax)
      {
        notused = 4;
        break;
      }
      if (si2 < 0)
        ++imag; /* commutator is odd */
      if (pvib1->wt[4] != pvib2->wt[4])
      {
        notused = 5;
        break;
      }
      lv1 = pvib1->lvqn;
      lv2 = pvib2->lvqn;
      if (iv12q == ilim && lv1 < 0)
        break;
      if (kavg > 0 && lv1 < 0)
      {
        if (lv1 != lv2)
          break;
        ++imag;
        isym = 3 - isym; /* 8x type operator */
      }
      if (ctx->glob.nofc && ifc != 0)
      {
        if (ifc < 0)
          ifc = -1 - ifc;
        kavg = ifc;
        ifc = 0;
        spar_now->fc = C0;
        spar_now->kavg = (char)kavg;
      }
      if (ODD(pvib1->lvupper))
      {
        --iv12;
        --iv1;
      }
      if (ODD(pvib2->lvupper))
      {
        iv12 -= ifac;
        --iv2;
      }
      /* setup for l-doubling */
      ldel = lv1 - lv2;
      if (ldel != 0)
      {
        if (ldel < 0)
          ldel = -ldel;
        if (ivdif < 0) /* DELTA K * DELTA L < 0 */
          ldel = -ldel;
        ctx->glob.lsym = FALSE;
      }
      else if (lv1 != 0 && (idflags & MLZ) == 0)
      {
        if (ctx->glob.newlz)
        {
          if (lv1 < 0)
          {
            /*  change symmetry for Lz operator */
            isym = 3 - isym;
            ++imag;
            idflags |= MLZ;
          }
        }
        else
        {
          if (lv1 < 0)
          {
            notused = 6;
            break;
          }
          i = (pvib1->lvupper ^ pvib2->lvupper) & 2;
          if (ODD2(isym + i))
          {
            /*  change symmetry for Lz operator */
            isym = 3 - isym;
            ++imag;
            idflags |= MLZ;
          }
        }
      }
      spar_now->mldel = (signed char)ldel;
      if (ODD(imag))
        idflags |= MODD;
      /* check for imaginary operators */
      if (iv1 == iv2 && ldel == 0 && isym != 0)
      {
        if (ODD(imag))
        {
          nimag[isym] += 1;
          idflags |= MIMAG;
        }
        else
        {
          nimag[isym + 4] += 1;
          idflags |= MIMAG;
        }
      }
      if (nitot < 3)
      { /* calculate alpha */
        if (gsym >= 3 && MOD(kd - ldel, gsym) != 0)
        {
          if (notused != 0)
            notused = 7;
          continue; /* matrix element breaks symmetry */
        }
        alpha = 0;
      }
      else
      {
        i = notused;
        if (notused != 0)
          notused = 8;
        alpha = MOD(ldel - kd, nitot);
        if (si1 <= ctx->itsym)
        {
          if (alpha != 0)
            continue;
          idflags |= MIDEN; /* identity spin operator */
        }
        if (alpha < 0)
          alpha += nitot;
        if (alpha == 0 || (alpha + alpha) == nitot)
        {
          if (ityi > 160)
            continue;
          if (ityi > 80 && kd == 0 && lv1 == 0 && lv2 == 0)
            continue;
        }
        else
        {
          ctx->glob.esym = FALSE;
          spar_now->zfac *= 0.5;
        }
        notused = i;
        if (ityi >= 80)
        {
          i = ityi / 80; /* i = 3,2,1,0 */
          if (ODD2(i))
          {
            if (ODD(isym + imag))
              spar_now->zfac = -spar_now->zfac;
            idflags |= MSYM2;
            isym = 3 - isym;
          }
          if (ODD(i))
          {
            alpha += nitot; /* special for quasi-diagonal */
          }
        }
      }
      if (ikq != 0 || neuler != 0 || sznz != 0 || ifc != 0 || kavg > 0)
        idflags |= MNOUNIT; /* no unit matrix */
      spar_now->alpha = (unsigned char)alpha;
      if (isym == 0 && iv1 == iv2 && kd != 0)
        iret = 1;
      if (testwt(ctx, pvib1, pvib2, isym, alpha))
      {
        if (notused != 0)
          notused = 9;
        continue;
      }
      if (ODD(neuler))
        isym = 3;
      ivsym = (unsigned int)isym + ((unsigned int)iv12 << 2);
      spar_now->ipsym = ivsym;
      spar_match = spar_last = NULL;
      for (spar_now = ctx->spar_head[iv2]; spar_now != NULL;
           spar_now = spar_now->next)
      {
        if (spar_now->ipsym <= ivsym)
          break;
        spar_last = spar_now;
      }
      if (EVEN(neuler))
      { /* not Euler denominator parameter */
        kl = idflags & MMASK;
        for (/* no init */; spar_now != NULL; spar_now = spar_now->next)
        {
          if (spar_now->ipsym < ivsym)
            break;
          spar_last = spar_now;
          klp = idpars(spar_now, &ikqp, &neulerp, &ltp, &ldp, &kdp, &insp,
                       &si1p, &si2p, &sznzp, &ifcp, &alphap, &ldelp, &kavgp);
          klp &= MMASK;
          /* check for same K dependence */
          if (ikq == ikqp && ld == ldp && kd == kdp && kl == klp &&
              sznz == sznzp && ifc == ifcp && alpha == alphap &&
              ldel == ldelp && neuler == 0 && neulerp == 0 && kavg == kavgp)
          {
            /*  check for same operator except power of N*(N+1) */
            if (lt == ltp && ins == insp && si1 == si1p && si2 == si2p)
            {
              idflags |= MNSQ;
            }
            else
            {
              if (TEST(idflags, MNSQ))
                break;
            }
            spar_match = spar_now;
          }
          else
          { /* operator with different K dependence */
            if (TEST(idflags, MNSQ))
              break;
          }
        } /* end loop over previous parameters */
      }
      notused = 0;
      if (spar_match != NULL)
      {
        if (ipar < 0)
          break; /* already have a unit operator */
        idflags |= MCOS_OK;
        spar_last = spar_match;
      }
      spar_now = spar_free;
      spar_free = NULL;
      if (spar_last != NULL)
      {
        spar_now->next = spar_last->next;
        spar_last->next = spar_now;
      }
      else
      {
        spar_now->next = ctx->spar_head[iv2];
        ctx->spar_head[iv2] = spar_now;
      }
      spar_last = NULL;
      spar_now->ip = ipar;
      spar_now->flags = (short)idflags;
    } while (ityi > 1);
    i = ipar;
    if (iv12q == ilim)
    {
      ++iv1d;
      if (iv1d >= ctx->glob.nvib)
        --ipar;
    }
    else
    {
      --ipar;
    }
    if (i > ipar)
    {
      if (notused > 0)
      {
        putbcd(ctx->sbcd, NSBCD, &idpar[ibcd]);
        if (notused * sizeof(char *) >= sizeof(strej))
          notused = 0;
        fprintf(lu,
                " WARNING: parameter %6d %s has no matrix elements: %s",
                (i + 1), ctx->sbcd, strej[notused]);
      }
      iv1d = 0;
      ibcd -= nbcd;
      notused = 0;
      ;
    }
  } /* end ipar loop */
  if (spar_free != NULL)
    free(spar_free);
  if (ctx->glob.esym)
    ctx->glob.lsym = FALSE;
  if (ctx->glob.stdphase != 0)
  { /* preset phase */
    ctx->glob.phasemask = 7;
    ctx->glob.stdphase &= 7;
  }
  else
  {
    ctx->glob.phasemask = 0;
    for (isym = 1; isym <= 3; ++isym)
    {
      /* check for phase change */
      if (nimag[isym] == 0 && nimag[isym + 4] == 0)
        continue;
      i = ctx->ipwr2[isym];
      ctx->glob.phasemask |= i;
      if (nimag[isym] > nimag[isym + 4])
      {
        ctx->glob.stdphase ^= i;
        i = nimag[isym];
        nimag[isym] = nimag[isym + 4];
        nimag[isym + 4] = i;
      }
    }
  }
  kk = ctx->glob.stdphase;
  if (kk > 0)
  {
    /* update phases */
    kk = kk ^ 5;
    for (isym = 1; isym <= 3; ++isym)
    {
      ctx->ixphase[isym] = kk & 1;
      kk = kk >> 1;
    }
  }
  nt = nimag[1] + nimag[2] + nimag[3];
  if (nt == 0)
    return iret;
  for (iv2 = 0; iv2 < ctx->glob.nvib; ++iv2)
  {
    spar_now = NULL;
    for (;;)
    {
      spar_last = spar_now;
      spar_now = (spar_now == NULL) ? ctx->spar_head[iv2] : spar_now->next;
      if (spar_now == NULL)
        break;
      isym = (int)(spar_now->ipsym & 7);
      if (isym == 0)
        continue;
      idflags = spar_now->flags;
      if ((idflags & MIMAG) == 0)
        continue;
      kk = 0;
      if (TEST(idflags, MODD))
        ++kk;
      if (TEST(ctx->glob.stdphase, ctx->ipwr2[isym]))
        ++kk;
      if (ODD(kk))
      { /* bad phase */
        ipar = spar_now->ip;
        ibcd = ipar * nbcd;
        putbcd(ctx->sbcd, NSBCD, &idpar[ibcd]);
        fprintf(lu,
                " WARNING: parameter %6d %s is imaginary and will not be used\n",
                (ipar + 1), ctx->sbcd);
        spar_free = spar_now->next;
        if (spar_last == NULL)
        {
          free(spar_now);
          ctx->spar_head[iv2] = spar_free;
        }
        else
        {
          free(spar_last->next);
          spar_last->next = spar_free;
        }
        spar_free = NULL;
        spar_now = spar_last;
        if (--nt <= 0)
          return iret;
      }
    }
  }
  return iret;
} /* pasort */

/**
 * @brief Gets the sub-block structure and key quantum numbers for a given F-block.
 * @param block_index Index of the F-block (0 to num_F_blocks-1).
 * @param f_qn_times_2_out Output: Pointer to store 2*F for this block.
 * @param stat_weights_block_out Output: Array to store statistical weights [A,E1,E2/B,gsym,nqn] for this block.
 * @param sub_block_pointers_out Output: Array to store start index of each sub-block (Wang*spin) within this F-block.
 * @param min_k_values_out Output: Array to store minimum K value for each sub-block.
 * @param vib_spin_sym_packed_out Output: Array to store packed (Vib,Sym,SpinPattern) identifier for each sub-block.
 * @return int Number of sub-blocks in this F-block. Returns 0 if block_index is invalid or no states.
 */
int getqq(struct SpinvContext *ctx, const int iblk, int *f, int *iwtb, short *sbkptr, short *kmin, int *vs)
{
  /* subroutine to get sub-block structure and major quantum numbers */
  /*     on entry: */
  /*         MXBLK=biggest block size */
  /*         IBLK= block number */
  /*     on return: */
  /*         GETQQ= number of sub-blocks */
  /*         F= F quantum number *2 */
  /*         IWTB= statistical weight */
  /*         SBKPTR= pointer to beginning index for sub-block */
  /*         KMIN= minimum K for each sub-block */
  /*         VS= vibration and symmetry for each sub-block */
  SVIB *pvinfo;
  short *jjs;
  int ibgn, iend, ksym, i, n, nsize, isblk, nsblk, k, kk, nn, iv;
  int mss, nset, nsym, ff;
  unsigned int mvs;

  isblk = iblk - 1;
  if (ctx->glob.nbkpj <= 0)
    return 0;
  ff = isblk / ctx->glob.nbkpj;
  isblk -= ctx->glob.nbkpj * ff; /* form remainder */
  ff = ff + ff;
  /*  get information for sub-blocks */
  nsize = nsblk = nsym = 0;
  ibgn = ctx->blkptr[isblk];
  iend = ctx->blkptr[isblk + 1];
  for (i = ibgn; i < iend; ++i)
  {
    mvs = (unsigned int)ctx->moldv[i];
    ksym = (int)mvs & 3;
    iv = (int)(mvs >> 2) & ctx->glob.msmask;
    mss = (int)(mvs >> ctx->glob.msshft);
    pvinfo = &ctx->vinfo[iv];
    if (i == ibgn)
      nsym = setgsym(ctx, (int)pvinfo->gsym);
    jjs = pvinfo->spt;
    nset = jjs[0];
    jjs += nset * mss;
    if (ff < jjs[nset - 1])
      continue;
    nn = jjs[0] + ff;
    n = nn >> 1;
    if (ODD(n))
      ksym = 3 - ksym;
    k = pvinfo->knmin[ksym];
    if (k < 0)
    {
      kk = MOD(jjs[ctx->itsym + 1] >> 2, nsym); /* get spin symmetry */
      k = ksym;
      if (ctx->is_esym[kk] < 0)
        k += 2;
      k &= 2;
    }
    kk = pvinfo->knmax;
    if (kk > n)
      kk = n;
    kk -= k;
    if (kk >= 0)
    {
      sbkptr[nsblk] = (short)nsize;
      kmin[nsblk] = (short)k;
      vs[nsblk] = (int)mvs;
      nsize += 1 + (kk >> 1);
      if (nsblk == 0)
      {
        *f = ff - (nn & 1);
        kk = -4;
        if (pvinfo->wt[4] != 0)
          kk -= jjs[ctx->itsym + 1];
        getwt(ctx, pvinfo, (int)(mvs & 3), kk, iwtb);
        iwtb[3] = pvinfo->gsym;
        iwtb[4] = pvinfo->nqn;
      }
      ++nsblk;
    }
  }
  if (nsblk > 0)
    sbkptr[nsblk] = (short)nsize;
  return nsblk;
} /* getqq */