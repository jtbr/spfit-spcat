/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 3 March 2004 */
/*   28 Sept.  2005: fix phase bug in getzitot for reduced matrix elements */
/* split from original spinv.c */

/* Utility functions, parsing, special calculations (intensity, direction cosine elements, n_k factors, quantum numbers) */
/*
idpari(): Parses the non-vibrational part of a BCD parameter ID.
idpars(): Unpacks a SPAR structure into individual field variables (less critical to move if only used by hamx).
ffcal(): Calculates N<sub>K</sub> factors for P<sup>±</sup> operators.
intens(): The main intensity calculation function.
dircos(): Calculates direction cosine matrix elements (crucial for both intensities and some Hamiltonian terms, but its primary role here is for transition moments).
getqs(): Gets quantum numbers for a sub-block basis state.
getqn(): Gets the full list of quantum numbers for a specific final state index.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "splib/calpgm_types.h"
#include "splib/blas_compat.h"
#include "splib/ulib.h"
#include "splib/cnjj.h"
#include "splib/slib.h"
#include "splib/catutil.h"
#include "spinit.h"
#include "spinv_internal.h"
#include "SpinvContext.hpp"
#include "common/CalError.hpp"

/**
 * @brief Parses the non-vibrational part of a BCD parameter identifier into a SPAR structure.
 *        It decodes fields for K^2 power, operator type, spin involvement, Fourier coefficients, etc.
 *        This function can be called multiple times for a single BCD ID if it represents
 *        multiple sub-operator types (e.g. for hybrid operators like P_a P_b + P_b P_a).
 * @param idval_bcd_param_id_no_vib Pointer to the BCD array for the non-vibrational part of the parameter ID. This array is NDECPAR long.
 * @param itp_subtype_in Input sub-parameter type. If 0, it's the first call for this BCD ID.
 *                       If >0, it's a subsequent call to get the next sub-operator component.
 * @param pspar_out Output: Pointer to the SPAR structure to be filled with parsed fields.
 * @return int The next itp_subtype_in value to use for further calls (if more sub-operators exist for this ID),
 *             or 0 if no more sub-operators for this BCD ID or if an invalid combination is encountered.
 */
int idpari(struct SpinvContext *ctx, bcd_t *idval, int itp, SPAR *pspar)
{
  /* subroutine to parse parameter identifier for complex parameter type */
  /*     on entry: */
  /*         IDVAL= parameter ID / vibrational field */
  /*         ITP= value of sub-parameter type for hybrid operator */
  /*     on return: */
  /*         pspar modified */

  /*0  1  2   3   4   5   6   7  8   9 */
  static int itpnxt[10] = {0, 0, 0, 2, 6, 6, 2, 0, 9, 0};
  static int itpop[10] = {0, 0, 0, 10, 40, 40, 10, 1, 0, 1};
  static int zfacv[10] = {1, 2, 6, 3, -8, 8, -6, 2, 4, -4};
  int ins, ity, nsx, isy, itysav;
  int itmp, i, iphaz, ldt, lnt, kdt, ksq, si1, si2, exval;
  bcd_t ibtmp;

  iphaz = 0;
  if (itp > 0)
  {
    /* get current sub-parameter type from previous */
    if (itp < 10)
    {
      itp = itpnxt[itp];
      if (itp == 0)
        return 0;
    }
    else
    {
      iphaz = itp >> 4;
      itp = itpnxt[itp & 15];
    }
  }
  pspar->zfac = 1.;
  pspar->euler = C0;
  pspar->msznz = C0;
  pspar->fc = C0;
  exval = 0;
  // The actual KSQ and NJQ are in idval[1] (which is idpar[ibcd + 2])
  // idval[0] contains the last two digits of the parameter ID
  ibtmp = idval[0]; // TODO: SHOULD read from the SECOND BCD element idval[1] (BUG in original code??!)
  nsx = (int)(ibtmp & 0x0f); // Lower nibble is NJQ
  ksq = (int)(ibtmp >> 4);   // Upper nibble is KSQ
  ity = bcd2i(idval[1]);
  itysav = ity * 10 + ksq;
  ibtmp = idval[2];
  ins = (int)(ibtmp & 0x0f);
  si1 = (int)(ibtmp >> 4);
  if (ins > 0 && ctx->nspin == 0)
    return 0;
  if (ins >= 5 && ctx->itsym == 0)
    return 0;
  ibtmp = idval[3];
  si2 = (int)(ibtmp & 0x0f);
  itmp = (int)(ibtmp >> 4);
  if (si2 > si1)
  {
    /*  make sure SI1 > SI2 */
    i = si1;
    si1 = si2;
    si2 = i;
    if (si2 == 0)
      si2 = -1;
  }
  if (si1 > ctx->nspin)
    return 0;
  if (ity > 90)
  {
    if (si1 != 0 || ins >= 5)
      return 0;
  }
  ibtmp = idval[4];
  if (ibtmp != (bcd_t)0)
  {
    if (ibtmp >= (bcd_t)0x60)
      return 0;
    itmp += 10 * (int)(ibtmp & 0x0f);
    exval = (int)ibtmp & 0xf0;
  }
  pspar->kavg = C0;
  if (itmp != 0)
  {
    i = (itmp - 1) / 10;
    itmp -= (i >> 1) * 10;
    if (ODD(i))
      itmp = 9 - itmp;
    pspar->fc = (char)itmp;
  }
  pspar->msi1 = (char)si1;
  pspar->msi2 = (char)si2;
  pspar->njq = (unsigned char)nsx;
  /*        exclude S * S and I * S from aa, bb, cc operators */
  if (itp == 2)
  {
    if (si2 > 0)
    {
      itp = itpnxt[itp];
      if (itp == 0 && iphaz == 0)
        return 0;
    }
  }
  if (itp != 0)
  {
    if (itp == 6 && si1 == 0)
    {
      itp = 8;
      pspar->njq = (unsigned char)(nsx + 1);
    }
  }
  else
  {
    if (iphaz != 0)
    {
      if (ins >= 5)
        iphaz -= 1;
      else if (ctx->glob.nitot >= 3)
        iphaz -= 5;
    }
    else
    {
      if (ins >= 5)
        iphaz = 4; /* iphaz = 4,3,2,1,0 */
      if (ctx->glob.nitot >= 3)
        iphaz += 15;
    }
    if (ity == 0)
    {
      itp = 1;
    }
    else if (ity <= 3)
    {
      itp = ity;
      if (ctx->glob.oblate)
        itp = ctx->revsym[itp];
      itp += 2;
      /*  ITP=3,4,5 for A,B,C */
      if (itp == 3 && si1 == 0)
        itp = 7;
    }
    else
    {
      itp = 1;
      itmp = ity;
      if (ity >= 92)
        itmp = (ity - 92) >> 1;
      /*  negate odd powers of P+**2 */
      if (itmp <= 11 && EVEN(itmp))
        pspar->zfac = -(pspar->zfac);
    }
  }
  isy = ity / 20;
  if (isy > 3)
  {
    isy = 0;
  }
  else if (isy > 0)
  {
    if (ctx->glob.oblate)
    {
      itysav += (ctx->revsym[isy] - isy) * 200;
    }
    else
    {
      isy = ctx->revsym[isy];
    }
  }
  pspar->ipsym = (unsigned int)isy;
  if (itp > 1)
  {
    itysav += itpop[itp] - ity * 10;
    pspar->zfac *= 2. / (double)zfacv[itp];
  }
  itp += iphaz << 4;
  if (ins >= 5)
  {
    ins -= 5;
    if (iphaz >= 5)
      iphaz = MOD(iphaz, 5);
    pspar->msznz = (unsigned char)(iphaz + 1);
  }
  pspar->mins = (unsigned char)ins;
  ity = itysav / 10;
  ksq = itysav - ity * 10;
  ldt = kdt = 0;
  if (ity <= 3)
  {
    if (ity != 0)
      ldt = 2;
  }
  else if (ity < 12)
  {
    ldt = kdt = ity + ity - 6;
  }
  else if (ity < 20)
  {
    kdt = ity + ity - 22;
    ldt = kdt + 1;
  }
  else if (ity < 80)
  {
    isy = ity / 20;
    ldt = kdt = ity - isy * 20 + 1;
    if (ODD((isy >> 1) + ldt))
      --kdt;
  }
  else if (ity < 90)
  {
    /* K energies coded as 8n  with K=10*n+KSQ */
    if (exval != 0)
      return 0;
    pspar->kavg = (char)(10 * (ity - 80) + ksq);
    ksq = 0;
  }
  else
  {
    if (exval != 0)
      return 0;
    i = ity - 90;
    if (ODD(i))
    {
      if (ksq > 1)
      { /* augmented Euler series */
        --i;
        ksq += 8;
        kdt = i;
      }
    }
    else
    {
      kdt = i;
    }
    pspar->euler = (unsigned char)(i + 2);
  }
  pspar->ksq = (char)ksq;
  if (si1 == 0)
  {
    lnt = 0;
  }
  else if (si2 <= 0)
  {
    if (isy == 0 || ldt == 0)
      si2 = 0;
    else if (si2 < 0)
      --ldt;
    lnt = 1;
  }
  else if (ldt < 2)
  {
    lnt = 0;
  }
  else
  {
    lnt = 2;
  }
  pspar->mln = (unsigned char)lnt;
  pspar->mkdel = (unsigned char)kdt;
  /*  try to use raising operator */
  if (ldt == kdt && ldt > lnt)
    ldt -= (ldt - lnt) & 0x7e; /* subtract even number */
  pspar->mld = (unsigned char)ldt;
  if (exval != 0)
  {
    /* check for Euler restrictions */
    ibtmp = (bcd_t)exval;
    exval = exval >> 3;
    if (idval[1] == (bcd_t)0 && idval[2] == (bcd_t)0 &&
        idval[3] == (bcd_t)0 && idval[4] == ibtmp)
    {
      /* energy Euler term */
      if (idval[0] == (bcd_t)0)
        return 0; /* bad energy-like term */
    }
    else if (exval == 2)
    {
      exval = 12; /* operator Euler term */
    }
    pspar->euler = (unsigned char)exval;
  }
  return itp;
} /* idpari */

/**
 * @brief Parses a parameter identifier (from SPAR struct) into its constituent fields.
 * @param pspar Pointer to the SPAR structure containing the packed parameter ID and pre-parsed fields.
 * @param ksq_out Output: power of K^2.
 * @param itp_out Output: Euler series type/index (pspar->euler).
 * @param tensor_order_N_out Output: tensor order in N. Renamed ln.
 * @param dir_cos_order_out Output: order of direction cosine part. Renamed ld.
 * @param k_delta_out Output: change in K. Renamed kdel.
 * @param n_dot_s_power_out Output: power of N.S. Renamed ins.
 * @param spin1_idx_out Output: index of first spin. Renamed si1.
 * @param spin2_idx_out Output: index of second spin. Renamed si2.
 * @param sznz_type_out Output: SzNz operator type. Renamed sznz.
 * @param fourier_coeff_idx_out Output: Fourier coefficient index. Renamed ifc.
 * @param alpha_sym_out Output: I_tot alpha symmetry component. Renamed alpha.
 * @param l_delta_out Output: change in l. Renamed ldel.
 * @param k_avg_out Output: K_avg value if applicable. Renamed kavg.
 * @param njq_out Output: Power of N(N+1) operator
 * @return int Flags from the pspar structure.
 */
int idpars(SPAR *pspar, int *ksq_out, int *itp_out,
           int *tensor_order_N_out, int *dir_cos_order_out, int *k_delta_out,
           int *n_dot_s_power_out, int *spin1_idx_out, int *spin2_idx_out,
           int *sznz_type_out, int *fourier_coeff_idx_out, int *alpha_sym_out,
           int *l_delta_out, int *k_avg_out, int *njq_out) /* Renamed parameters for clarity */
{
  /*     subroutine to parse parameter identifiers to subfields */ /* Original comment */
  /*     on entry: */ /* Original comment */
  /*         IDPAR= 10*( parameter ID/1000) + ITP */ /* Original comment - this refers to an older BCD packing, now fields are directly in pspar */
  /*     on return: */ /* Original comment */
  /*         KSQ= power of K*K */ /* Original comment */
  /*         ITP= parameter subtype */ /* Original comment - now pspar->euler is primary type, itp_out is that */
  /*         LN= tensor order */ /* Original comment */
  /*         LD= order of direction cosine ( L <= LD ) */ /* Original comment */
  /*         KDEL= change in K */ /* Original comment */
  /*         INS=power of N.S */ /* Original comment */
  /*         SI1,SI2=spin identifiers */ /* Original comment */
  int flags_return; /* Renamed iret */
  *ksq_out = (int) pspar->ksq;
  *itp_out = (int) pspar->euler; /* This is the primary Euler type/flag now */
  *tensor_order_N_out = (int) pspar->mln;
  *dir_cos_order_out = (int) pspar->mld;
  *k_delta_out = (int) pspar->mkdel;
  *n_dot_s_power_out = (int) pspar->mins;
  *spin1_idx_out = (int) pspar->msi1;
  *spin2_idx_out = (int) pspar->msi2;
  *sznz_type_out = (int) pspar->msznz;
  *fourier_coeff_idx_out = (int) pspar->fc;
  *alpha_sym_out = (int)pspar->alpha;
  *l_delta_out = (int)pspar->mldel;
  *k_avg_out = (int)pspar->kavg;
  *njq_out = (int)pspar->njq;
  flags_return = (int)pspar->flags;
  return flags_return;
} /* idpars */

/**
 * @brief Calculates factors related to N_K = sqrt(N(N+1)-K(K-1)) for P+- operators.
 *        Stores sqrt( N(N+1) - K(K-1) ) for K from 0 up to k_val_limit.
 *        Used by dircos when operator involves P+ or P- type terms (DK_op > L_op).
 * @param n_val_for_factor N quantum number.
 * @param k_val_limit Maximum K value for which to calculate the factor.
 * @param output_factors_array Output: Array to store the calculated factors ff[K] = sqrt(N(N+1)-K(K-1)).
 * @return int Always 0.
 */
int ffcal(const int nff, const int kff, double *ff)
{
  static int nlast = -1;
  static int klast = 0;
  static double sq;
  int k;

  if (nff != nlast)
  {
    nlast = nff;
    sq = (double)nff;
    sq *= nff + 1;
    k = 0;
  }
  else
  {
    if (klast >= kff)
      return 0;
    k = klast;
  }
  klast = kff;
  if (klast > nff)
    klast = nff;
  while (k < klast)
  {
    ff[k] = sqrt(sq);
    ++k;
    sq -= k + k;
  }
  return 0;
} /* ffcal */

/**
 * @brief Calculates intensity matrix elements for transitions between states in different blocks.
 * @param iblk_upper_state Upper state F-block index.
 * @param isiz_upper_dim Dimension of upper state F-block.
 * @param jblk_lower_state Lower state F-block index.
 * @param jsiz_lower_dim Dimension of lower state F-block.
 * @param num_dipole_params Number of dipole moment parameters.
 * @param bcd_dipole_ids Array of dipole moment identifiers in BCD format.
 * @param dipole_values Array of dipole moment values.
 * @param intensity_matrix_out Output intensity matrix S_ij = <bra_i|mu|ket_j>.
 * @return int (2*F_upper + 1) if transitions are allowed and calculated, 0 otherwise.
 *         Effectively, the degeneracy of the upper F state.
 */
int intens(struct SpinvContext *ctx, const int iblk, const int isiz, const int jblk,
           const int jsiz, const int ndip, const bcd_t *idip, const double *dip, double *s)
{
  /*  function to calculate intensities */
  /*     on entry: */
  /*         IBLK,JBLK= block numbers */
  /*         ISIZ,JSIZ= block sizes */
  /*         IDIP= dipole id codes , IDIP(0)= number of dipoles */
  /*         DIP= dipole values */
  /*     on return: */
  /*         S= transition dipole matrix in original basis */
  /*         INTENS returns true is S is not 0 */
  static int idipoff;
  SDIP *pdip;
  double dd;
  long ndms;
  unsigned int ijv;
  int ncos, ifup, i, k, lv, n, icase, ibase, mkd, ibcd, nbcd, isunit;
  int jbase, kbgni, kbgnj, nblki, nblkj, ixtmp, iret, kd, ld, npair, ldel;
  int ix, jx, nx, si1, si2, iff, jff, nni, nnj, iwt[5], jwt[5], ndecv;
  int ivv, jvv, ixx, jxx, ifc, ksym, kl, ioff, joff, alpha, nitot, dipoff;
  bcd_t bijv1, bijv2, bijv3;

  if (ctx->ndmx <= 0)
    throw NumericError("working vectors not allocated", CalErrorCode::WorkingVectorTooShort);
  dipoff = idipoff;
  if (dipoff >= ctx->nddip)
    dipoff = 0;
  idipoff = dipoff + ndip;
  iret = 0;
  ctx->cgetv[0].cblk = 0;
  ctx->cgetv[1].cblk = 0;
  /*     get quantum information */
  nblki = getqq(ctx, iblk, &iff, iwt, ctx->ibkptr, ctx->ikmin, ctx->ivs);
  ioff = ctx->glob.maxblk;
  nblkj = getqq(ctx, jblk, &jff, jwt, &(ctx->ibkptr)[ioff], &(ctx->ikmin)[ioff], &(ctx->ivs)[ioff]);
  ixx = iff - jff;
  if (ixx < -2 || ixx > 2)
    return 0;
  ixx = iff + jff;
  if (ixx < 2 || ODD(ixx))
    return 0;
  if (iwt[3] != jwt[3])
    return 0;
  if (checkwt(ctx, iwt, jwt) != 0)
    return 0;
  nitot = ctx->glob.nitot;
  alpha = 0;
  bijv3 = (bcd_t)0;
  /* clear dipole matrix */
  dclr(isiz, jsiz, s, 1);
  /* loop over sub-blocks */
  ndms = isiz;
  ndecv = ctx->glob.vibdec;
  nbcd = (int)idip[0] & 0x7f;
  for (ixx = 0; ixx < nblki; ++ixx)
  {
    ibase = ctx->ibkptr[ixx];
    n = ctx->ibkptr[ixx + 1] - ibase;
    kbgni = ctx->ikmin[ixx];
    ixtmp = ctx->ivs[ixx];
    getqs(ctx, ixtmp, iff, n, kbgni, ctx->ixcom, ctx->iscom, &ivv);
    nni = ctx->iscom[0];
    for (jxx = 0; jxx < nblkj; ++jxx)
    {
      joff = jxx + ioff;
      jbase = ctx->ibkptr[joff];
      n = ctx->ibkptr[joff + 1] - jbase;
      kbgnj = ctx->ikmin[joff];
      ixtmp = ctx->ivs[joff];
      getqs(ctx, ixtmp, jff, n, kbgnj, ctx->jxcom, ctx->jscom, &jvv);
      nnj = ctx->jscom[0];
      nx = ((nni > nnj) ? nni : nnj) >> 1;
      i = (ivv < jvv) ? ivv : jvv;
      ijv = (unsigned int)(ivv + jvv + i * ctx->glob.vibfac);
      ijv = (ijv << 2) + blksym(ctx->ixcom, ctx->jxcom);
      bijv1 = (bcd_t)ijv;
      bijv2 = (bcd_t)(ijv >> 8);
      if (ndecv == 3)
        bijv3 = (bcd_t)(ijv >> 16);
      for (i = 0, ibcd = 0; i < ndip; ++i, ibcd += nbcd)
      {
        /*  set dipole for transition */
        if (bijv1 != idip[ibcd + 1])
          continue;
        if (bijv2 != idip[ibcd + 2])
          continue;
        if (bijv3 != idip[ibcd + 3])
          continue;
        pdip = &(ctx->dipinfo)[dipoff + i];
        kl = pdip->flg;
        if (TEST(kl, MINOQ) && nni == nnj)
          continue;
        /* isym = (int)(bijv1 & 3); */
        si1 = (int)(idip[ibcd + 4] >> 4) & 0x0f;
        si2 = 0;
        icase = bcd2i(idip[ibcd + 5]);
        kd = (int)pdip->kd;
        ifc = (int)pdip->fc;
        ldel = (int)pdip->ldel;
        lv = 1;
        ld = (int)pdip->ld;
        if (si1 > 0)
          lv = TEST(kl, MELEC) ? 2 : 0;
        ksym = 0;
        if (nitot >= 3 && kd > 0)
        {
          alpha = nitot;
          ksym = 1;
        }
        do
        {
          if (ksym == 0)
          {
            alpha = 0;
          }
          /*  find matrix elements */
          mkd = getmask(ctx, ctx->ixcom, ctx->jxcom, kd, ldel, kl, alpha);
          if (mkd == 0)
            continue;
          npair = getll(ctx, 2, ld, lv, kd, si1, si2, ctx->lscom, ctx->iscom, ctx->jscom);
          if (npair < 0)
            continue;
          ifup = TEST(kl, MDIPI) ? -1 : 1;
          isunit = (int)pdip->kavg;
          if (isunit > nx)
            continue;
          ncos = dircos(ctx, ctx->ixcom, ctx->jxcom, ld, kd, ctx->ndmx, ctx->wk, ctx->idx, ctx->jdx,
                        ifup, kl, mkd, &isunit);
          if (ncos == 0)
            continue;
          if (ncos < 0)
            throw NumericError("DIRCOS WORKING VECTOR TOO SHORT IN INTENS", CalErrorCode::WorkingVectorTooShort);
          dd = dip[i] * pdip->fac;
          /* correct for special cases */
          switch (icase)
          {
          case 3:
            if (nni == nnj)
            {
              k = nx * (nx + 1);
            }
            else
            {
              k = nx * nx;
            }
            dd *= k;
            break;
          case 4:
            symksq(1, kbgni, kbgnj, ncos, ctx->wk, ctx->idx, ctx->jdx);
            break;
          case 5:
          case 11:
            if (ifup > 0)
              break;
            if (nni < nnj)
              dd = -dd;
            dd *= nx;
            break;
          case 9:
            if (nni == nnj)
            {
              k = nx * (nx + 1);
              dd *= k;
            }
            else
            {
              k = nx * nx;
              dd *= k + 1;
            }
            dd *= k;
            break;
          case 10:
            dd *= nx * nx;
            break;
          case 7:
          case 8:
          case 12:
            specfc(ctx, ifc, ivv, jvv, kd, kbgni, kbgnj, ncos, ctx->wk, ctx->idx, ctx->jdx);
            break;
          }
          /*  correct for reduced matrix of N */
          if (ld != lv)
            dd *= rmatrx(ld, lv, ctx->ixcom, ctx->jxcom);
          /*  couple dipoles through the spins */
          tensor(ctx, &dd, ctx->iscom, ctx->jscom, ctx->lscom, ctx->ismap, npair, alpha);
          for (k = 0; k < ncos; ++k)
          {
            ix = ctx->idx[k] + ibase;
            jx = ctx->jdx[k] + jbase;
            s[ix + jx * ndms] += dd * ctx->wk[k];
          }
          iret = iff + 1;
        } while (--ksym >= 0);
      }
    }
  }
  return iret;
} /* intens */

/**
 * @brief Calculates direction cosine matrix elements <bra|Phi_Fgz|ket> in a K-basis.
 *        Handles various tensor orders (L=0,1,2) and K changes (DeltaK).
 *        Also incorporates l-doubling and I_tot symmetry effects via operator_flags and selection_mask.
 * @param xbra_qns Quantum numbers for the bra state.
 * @param xket_qns Quantum numbers for the ket state.
 * @param dir_cos_order Tensor order L of the direction cosine operator.
 * @param k_delta_op Overall Delta K for the operator.
 * @param max_elements_out Max number of elements that can be stored in output arrays.
 * @param dir_cos_elements_out Output: Array to store the calculated direction cosine matrix elements (K-dependent part).
 * @param ibra_indices_out Output: Array of bra K-indices (relative to K_begin_bra) for each element.
 * @param iket_indices_out Output: Array of ket K-indices (relative to K_begin_ket) for each element.
 * @param ifup_triangle_flag Controls calculation for upper/lower triangle. If <0 (electric dipole), affects phase.
 * @param operator_flags Bit flags for operator properties (MODD, MLZ from SPAR/SDIP).
 * @param selection_mask 3-bit mask from getmask() indicating allowed K connections (K'=K+DK, K'=DK-K, K'=K-DK).
 * @param is_unit_matrix_type_ptr Output/Input: If input k_avg_val > 0, it's a hint for unit matrix.
 *                                Output is 1 if operator is effectively a unit matrix (K-independent), 0 otherwise.
 * @return int Number of non-zero direction cosine elements calculated. Returns negative of deficit if max_elements_out is too small.
 */
int dircos(struct SpinvContext *ctx, const int *xbra, const int *xket, const int ld, const int kd,
           const int ncmax, double *direl, short *ibra, short *iket, const int ifup, const int loff,
           int mask, int *isunit)
{ /*  SUBROUTINE TO  CALCULATE DIRECTION COSINE ELEMENTS */
/* XBRA,XKET INTEGER VECTORS CONTAINING: */
/*           0: DIMENSION OF SUB-BLOCK */
/*           1: SYMMETRY CODE FOR BLOCK (0=A,1=Bx,2=By,3=Bz) */
/*           2: 'N' QUANTUM NUMBER */
/*           3: BEGINNING K */
/*           4: L VALUE */
/*           5: VIBRATIONAL QUANTA */
/*           6: SPIN SYMMETRY */
/* L IS TENSOR ORDER OF COSINE */
/* KD IS TOTAL K QUANTIUM CHANGE */
/* DIREL IS A VECTOR OF DIRECTION COSINE ELEMENTS */
/* NCMAX IS MAXIMUM LENGTH OF DIREL */
/* IBRA,IKET INDICES OF MATRIX ELEMENTS IN DIREL (START WITH ZERO) */
/* IFUP =0 IF NO UPPER TRIANGLE TO BE CALC. */
/* (LOFF AND MSYM2) IF OPERATOR HAS ALTERNATE SYMMETRY UNDER ITOT */
/* (LOFF AND MLZ) IF OPERATOR INCLUDES LZ */
/* (LOFF AND MODD) IF OPERATOR IS ODD-ORDER */
/*  L=0  KD=0    DIREL=1. */
/* IFUP <0 IF ELECTRIC DIPOLE */
/* ODD ORDER PARAMETERS ARE IMAGINARY AND */
/* EVEN ARE REAL IF NOT ELECTRIC DIPOLE */
/* FOR ELECTRIC DIPOLE ODD ORDER IS REAL AND EVEN ARE IMAGINARY  */
/*  L=1  KD=0    DIREL=PHI(Z) */
/*  L=1  KD=1    DIREL=PHI(X) OR PHI(Y) */
/*  L=2  KD=0    DIREL=0.5*(3* PHI(Z)**2 -1.) */
/*  L=2  KD=1    DIREL=PHI(Z)*PHI(X)+PHI(X)*PHI(Z)  OR */
/*                      =PHI(Z)*PHI(Y)+PHI(Y)*PHI(Z) */
/*  L=2  KD=2    DIREL=2*(PHI(X)**2 - PHI(Y)**2) OR */
/*                      =PHI(X)*PHI(Y)+PHI(Y)*PHI(X) */
/*     IF KD > L THEN (P-)**(KD-L) OPERATOR MULTIPLIED BY */
/*               L,L OPERATOR */
/*     WANG FUNCTIONS MULTIPLIED BY i**(SYMMETRY) TO MAKE ALL */
/*               ODD ANGULAR MOMENTA IMAG. AND ALL EVEN ONES REAL */
/* IFC = FOURIER SERIES, > 0 IS COS, < 0 IS SIN */
/* MASK = 3-BIT MASK TO ENABLE MATRIX ELEMENT SUB-TYPE */
/*     BEST RESULTS ON SPEED IF N FOR BRA CHANGES SLOWEST */
/***********************************************************************/
#define NFAC_DIRCOS 20
  static double ff[MAXN_DIRCOS], fac[NFAC_DIRCOS];
  double ele, elex, eletmp;
  int kdel, kbra, kbgn, nbgn, kket, isgn, ncase0;
  int ncase, isum, kbra0, kdel2, kket0, kket2, i, k;
  int n, ikdel, kkdel, isbra, nnbra, nsbra, isket, idif;
  int ii, il, kk, lld, ir, ks, nnket, nsket, isrev, nval;
  int iphaz, isym, nbra, nket, kbra0x, kket0x, kavg, kbit;
  BOOL matsym;

  kbra0 = xbra[XKBGN];
  kbra0x = xbra[XIQN] & 1;
  kket0 = xket[XKBGN];
  kket0x = xket[XIQN] & 1;
  nsbra = xbra[XDIM];
  isbra = xbra[XSYM];
  nsket = xket[XDIM];
  isket = xket[XSYM];
  isym = (isbra ^ isket) & 3;
  if ((loff & MODD) != 0)
  {
    if (ncmax == 0)
    { /* initialize */
      fac[0] = 0.5;
      for (lld = 1; lld < NFAC_DIRCOS; ++lld)
      {
        ele = lld - 0.5;
        fac[lld] = fac[lld - 1] * sqrt(lld / ele);
      }
      fac[0] = 1. / fac[1];
      return 0;
    };
    isym = 3 - isym;
  }
  ncase = 0;
  kavg = (*isunit);
  *isunit = 0;
  /****** set up state phases, isum is power of i */
  isrev = isym >> 1;
  /* isym is odd if operator is imaginary */
  isum = ctx->ixphase[isket] - ctx->ixphase[isbra] + (isym & 1);
  if (ifup < 0) /* correct to make dipoles real */
    ++isum;
  if (ODD(isum) && xbra[XVIB] < xket[XVIB])
    isum += 2; /* imaginary vibration phase assumed */
  isgn = isum >> 1;
  if ((loff & MLZ) != 0 && xket[XLVAL] < 0)
    ++isgn; /* Lz operator found */
  /*********** TRY SCALAR ******************************/
  if (kd == 0 && ld == 0 && (mask & 5) != 0)
  {
    mask &= 2;
    ii = 0;
    ncase = nsket;
    n = nsbra;
    kdel = kbra0 - kket0;
    kk = kbra0;
    if (kdel != 0)
    {
      if (kdel > 0)
      {
        kdel = kdel >> 1;
        ncase -= kdel;
        if (ncase < 0)
          ncase = 0;
        kk = kbra0;
      }
      else
      {
        ii = (-kdel) >> 1;
        kdel = -ii;
        n -= ii;
        if (n < 0)
          n = 0;
        kk = kket0;
      }
    }
    if (ncase > n)
      ncase = n;
    if (ncase > 0)
    {
      idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      if (kavg > 0)
      { /* single K value */
        k = kavg - kk;
        if (k < 0 || ODD(k))
          return 0;
        k = k >> 1;
        if (k >= ncase)
          return 0;
        ii += k;
        kk += k;
        ncase = 1;
      }
      if (kk != 0)
      { /* kk = max beginning k */
        mask = 0;
      }
      else if (kbra0x == kket0x)
      {
        if (kbra0x == 0)
          mask = 0;
        kk = 2;
      }
      else
      { /* kk = 0 */
        mask = 0;
      }
      ele = 1.;
      if (ODD(isgn + isrev))
        ele = -ele;
      *isunit = 1;
      for (k = 0; k < ncase; ++k)
      {
        direl[k] = ele;
        ibra[k] = (short)ii;
        iket[k] = (short)(ii + kdel);
        ++ii;
      }
      if (kk == 0)
      {
        direl[0] *= fac[0];
        *isunit = 0;
      }
    }
    if (mask == 0)
      return ncase; /* no quasi-diagonal delta K = 0 */
  } /************** end Scalar *********/
  nbra = xbra[XNVAL];
  nket = xket[XNVAL];
  kdel = kd;
  ikdel = kdel - ld;
  if (ikdel <= 0)
  {
    kkdel = kdel;
  }
  else
  {
    k = kbra0 + ikdel + ((nsbra - 1) << 1);
    ffcal(nbra, k, ff);
    kkdel = ld;
  }
  kdel2 = kkdel + kkdel;
  lld = ld + ld;
  nnbra = nbra + nbra;
  nnket = nket + nket;
  if (ld > 0)
  {
    ele = (double)(nnbra + 1);
    if (nnket != nnbra)
      ele = sqrt(ele * (nnket + 1));
    if (kkdel != 0)
      ele *= fac[ld];
    isgn += nbra + kket0 + kkdel;
  }
  else
  {
    ele = 1.;
    if (kdel != 0)
      ele = 0.5;
  }
  matsym = TRUE;
  k = mask & 5;
  if (kdel == 0)
  {
    if (k != 0)
      mask = (mask & 2) + 1;
    matsym = FALSE;
  }
  else
  {
    if (k == 1 || k == 4)
      matsym = FALSE;
  }
  if (matsym)
  {
    for (k = 1; k <= 5; ++k)
    {
      if (xbra[k] != xket[k])
      {
        matsym = FALSE;
        break;
      }
    }
  }
  if (ODD(isgn + isrev))
    ele = -ele;
  kbit = kavg + kavg - (kdel & 1);
  ncase0 = 0;
  for (iphaz = 1; iphaz <= 4; iphaz = iphaz << 1)
  {
    /* iphaz = 1,2,4 */
    if (mask < iphaz)
      break;
    if (iphaz != 2)
    {
      /* iphaz = 1: +KDEL, iphaz = 4: -KDEL */
      /*  START CALCULATION OF  <K ....K-KDEL> */
      /*     KBRA=KBRA0+ (IL+ i)*KINC   IL >= 0 */
      /*     KKET=KKET0+ (IR+ i)*KINC   IR >= 0, i=0,... */
      /*     KBRA=KKET+KDEL */
      if ((mask & iphaz) == 0)
        continue;
      if (iphaz == 4)
      {
        if (matsym)
          break;
        kdel = -kdel;
        kkdel = -kkdel;
      }
      isum = kdel + kket0 - kbra0; /* ISUM = (IL-IR)*KINC */
      if (isum >= 0)
      {
        il = isum >> 1;
        ir = 0;
        kbra = kdel + kket0;
        kket = kket0;
      }
      else
      {
        il = 0;
        ir = (-isum) >> 1;
        kbra = kbra0;
        kket = kbra0 - kdel;
      }
      n = nsket - ir;
      i = nsbra - il;
      if (n > i)
        n = i;
      if (n <= 0)
        continue;
      if (kavg > 0)
      {
        kk = (kbit - kbra - kket) >> 1;
        if (kk < 0 || ODD(kk))
          continue;
        kket += kk;
        kbra += kk;
        kk = kk >> 1;
        if (kk >= n)
          continue;
        il += kk;
        ir += kk;
        n = 1;
      }
      nbgn = ncase;
      ncase += n;
      idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      kket2 = kket + kket;
      ncase0 = n;
      if (ld == 0)
      {
        for (n = nbgn; n < ncase; ++n)
        {
          direl[n] = ele;
          ibra[n] = (short)il;
          iket[n] = (short)ir;
          ++il;
          ++ir;
        }
      }
      else
      {
        for (n = nbgn; n < ncase; ++n)
        {
          kk = -kdel2 - kket2;
          direl[n] = ele * c3jj(nnbra, lld, nnket, kk, kdel2, kket2);
          kket2 += 4;
          ibra[n] = (short)il;
          iket[n] = (short)ir;
          ++il;
          ++ir;
        }
      }
      il = ir = 1;
      if (kbra == 0)
        il = kbra0x;
      if (kket == 0)
        ir = kket0x;
      if (ODD(il + ir)) /*  CORRECT FOR K=0 */
        direl[nbgn] *= fac[0];
      if (ikdel > 0)
      { /* ADD ON P- OR P+ OPERATOR TO LEFT */
        kbgn = kbra;
        if (kdel > 0)
          kbgn -= ikdel;
        for (n = nbgn; n < ncase; ++n)
        {
          eletmp = direl[n];
          kk = kbgn + ikdel;
          for (k = kbgn; k < kk; ++k)
          {
            eletmp *= ff[k];
          }
          direl[n] = eletmp;
          kbgn += 2;
        }
      }
    }
    else
    { /* iphaz == 2 */
      /*  DO QUASI-DIAGONAL PART */
      /*   <KDEL-K ..... K>   K=1,...,KDEL-1 */
      /*       KBRA=KBRA0+ (IL - i)*KINC  IL <= n-1 */
      /*       KKET=KKET0+ (IR + i)*KINC  IR >= 0  , i=0 .. n-1 */
      /*       KDEL=KKET+KBRA */
      kdel2 = -kdel2;
      if (ODD(isrev))
        ele = -ele;
      if ((mask & 2) == 0)
        continue;
      isum = kdel - kbra0 - kket0;
      if (isum < 0)
        continue;
      isum = isum >> 1; /* isum = (IL + IR) */
      ir = 0;           /* IR = 0 */
      if (isum >= nsbra)
        ir = isum - nsbra + 1;
      else if (kket0x == 0)
        ++ir;
      n = isum + 1;
      if (n > nsket)
        n = nsket;
      else if (kbra0x == 0)
        --n;
      n -= ir;
      if (n <= 0)
        continue;
      il = isum - ir;
      if (ifup == 0)
        n = (il - ir + 2) >> 1;
      kket = kket0 + (ir << 1);
      kbra = kdel - kket;
      nbgn = ncase;
      ncase += n;
      idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      if (ODD2(isbra + nnbra))
      {
        elex = -ele;
      }
      else
      {
        elex = ele;
      }
      n = nbgn;
      isgn = loff & MFCM;
      for (nval = nbgn; nval < ncase; ++nval)
      {
        kk = 0;
        if (kavg > 0)
        {
          kk = kket - kbra;
          if (kk < 0)
            kk = -kk;
          kk -= kbit;
        }
        else if (kbra == kket)
        {
          kk = isgn; /* skip if sin(0) */
        }
        if (kk == 0)
        {
          eletmp = elex;
          if (isgn != 0 && kbra > kket)
            eletmp = -eletmp;
          if (ld != 0)
          {
            kket2 = kket + kket;
            ks = -kdel2 - kket2;
            eletmp *= c3jj(nnbra, lld, nnket, ks, kdel2, kket2);
          }
          for (kk = 0; kk < ikdel; ++kk)
          {
            k = kbra - kk;
            if (k <= 0)
            {
              eletmp *= ff[-k];
            }
            else
            {
              eletmp *= ff[k - 1];
            }
          }
          direl[n] = eletmp;
          ibra[n] = (short)il;
          iket[n] = (short)ir;
          ++n;
        }
        --il;
        ++ir;
        kbra -= 2;
        kket += 2;
      } /* loop over elements */
      ncase = n;
    }
  } /* loop over 3 cases */
  if (matsym && ncase0 > 0 && ifup != 0)
  {
    k = ncase;
    ncase += ncase0;
    idif = ncmax - ncase;
    if (idif < 0)
      return idif;
    for (i = 0; i < ncase0; ++i, ++k)
    {
      direl[k] = direl[i];
      ibra[k] = iket[i];
      iket[k] = ibra[i];
    }
  }
  return ncase;
} /* dircos */

/**
 * @brief Gets quantum numbers for a given sub-block basis state.
 * @param packed_vib_spin_sym Packed (Vibrational_index, Symmetry, Spin_pattern_index). Renamed im.
 * @param f_qn_times_2 2*F quantum number for the block. Renamed iff.
 * @param sub_block_dim Dimension of this specific sub-block. Renamed nsiz.
 * @param k_begin Starting K value for this sub-block. Renamed kbgn.
 * @param ixcom_out Output array (NDXCOM size) to store rotational quanta. Renamed ixcom.
 * @param iscom_out Output array (MAXNS size) to store coupled angular momenta (2N, 2J, 2F1...). Renamed iscom.
 * @param vib_idx_out Output vibrational state index. Renamed iv.
 * @return int Spin pattern index (ispx) for this state.
 */
int getqs(struct SpinvContext *ctx, const int im, const int iff, const int nsiz,
          const int kbgn, int *ixcom, int *iscom, int *iv)
{
  /*     subroutine to get quanta corresponding to sub-block */
  /*     on entry: */
  /*         IM= position in sub-block mold */
  /*         IFF = F quantum number *2 */
  /*         NSIZ= size of sub-block */
  /*         KBGN= first K value */
  /*     on return: */
  /*         IXCOM = list of rotational containing */
  /*                 0: DIMENSION OF SUB-BLOCK */
  /*                 1: SYMMETRY CODE FOR BLOCK (0=A,1=Bx,2=By,3=Bz) */
  /*                 2: 'N' QUANTUM NUMBER */
  /*                 3: BEGINNING K */
  /*                 4: L VALUE */
  /*                 5: VIBRATIONAL QUANTA */
  /*                 5: ITOT SYMMETRY */
  /*                 5: EXTRA ITOT QUANTA */
  /*         ISCOM = list of angular momenta */
  /*         IV = vibrational quantum number   if L=0,1 */
  /*         IV = vibrational quantum number-1 if L=-1 */
  /*     find list of angular momenta */
  SVIB *pvinfo;
  int *isscom;
  short *iis, *jjs;
  int i, iiv, nset, nspinv, ispx, k, lupper;
  unsigned int mvs;

  mvs = (unsigned int)im;
  ispx = (int)(mvs >> ctx->glob.msshft);
  iiv = (int)(mvs >> 2) & ctx->glob.msmask;
  pvinfo = &(ctx->vinfo)[iiv];
  iis = pvinfo->spt;
  nset = iis[0];
  jjs = &iis[nset * ispx];
  i = (int)(mvs & 3);
  iscom[0] = iff + (*jjs);
  ixcom[XNVAL] = iscom[0] >> 1; /* find N */
  ixcom[XDIM] = nsiz;
  ixcom[XSYM] = i; /* find symmetry */
  if (nsiz < 0)
    setgsym(ctx, (int)pvinfo->gsym);
  k = kbgn;
  if (k < 0)
  {
    k = pvinfo->knmin[i];
    i = pvinfo->knmin[3 - i];
    if (k > i)
      k = i;
    if (k < 0)
      k = 0;
  }
  ixcom[XKBGN] = k;
  ixcom[XVIB] = iiv;
  ixcom[XISYM] = 0;
  lupper = pvinfo->lvupper;
  if (ODD(lupper))
    --iiv;
  *iv = iiv;
  i = pvinfo->lvqn;
  ixcom[XLVAL] = i;
  if (i != 0 || k != 0)
    i = 1;
  ixcom[XIQN] = i;
  if (ctx->nspin == 0)
    return 1;
  nspinv = nset - 1;
  if (nspinv > ctx->nspin)
    nspinv = ctx->nspin;
  isscom = &iscom[ctx->nspin];
  isscom[0] = iscom[0];
  for (i = 1; i < nspinv; ++i)
  {
    iscom[i - 1] = iff + jjs[i];
    isscom[i] = iis[i];
  }
  iscom[i - 1] = iff;
  isscom[i] = iis[i];
  while (i < ctx->nspin)
  {
    ++i;
    iscom[i - 1] = iff;
    isscom[i] = 0;
  }
  if (ctx->glob.nitot > 0)
  {
    if (ctx->glob.nitot >= 3)
    {
      if (ctx->itsym < nspinv)
      {
        i = jjs[ctx->itsym + 1];
        iscom[ctx->itsym] = i;
        k = MOD(i >> 2, ctx->glob.nitot);
        i -= k << 2;
        if (ctx->is_esym[k] != 0)
          ++i;
        ixcom[XISYM] = k;
        ixcom[XIQN] |= i;
      }
      else
      {
        iscom[ctx->itsym] = 0;
      }
    }
    iscom[ctx->itptr] = jjs[ctx->itptr + 1];
  }
  return ispx;
} /* getqs */

/**
 * @brief Gets the full set of quantum numbers for a specific state within a block.
 * @param f_block_idx F-block index. Renamed iblk.
 * @param state_idx_in_block Index of the state within the F-block (0 to Nblock-1). Renamed indx.
 * @param max_qn_fields Max number of quantum number fields to fill in qn_list_out. Renamed maxqn.
 * @param qn_list_out Output array to store the quantum numbers. Renamed iqn.
 * @param degeneracy_out Output: Pointer to store the degeneracy of the state. Renamed idgn.
 * @return int Returns a code related to Wang block sorting or K-parity (ncod).
 *             Positive for one K-parity, negative for the other, for asymmetric tops.
 *             1 for symmetric top / linear K=0, or if K-parity not distinguished.
 */
int getqn(struct SpinvContext *ctx, const int iblk, const int indx, const int maxqn, short *iqn, int *idgn)
{
  /*   subroutine to get quantum numbers */
  /*     on entry: */
  /*         IBLK= block number */
  /*         INDX= index within block */
  /*     on return: */
  /*         IQN= list of quantum numbers */
  /*         IDGN= degeneracy */
  /*         NCOD= seraching aid to find state within a Wang block */
  static int last;
  int ibgn, iend, iqsp, isym, i, k, nsblk, ioff, lower, iv, iqf;
  int n, ncod, ldgn, ix, kq, ivbase, lv;
  short mgsym, ngsym;

  /*  get sub-block info from store or GETQQ */
  /*     look at last used then oldest */
  if (iblk == ctx->cgetv[0].cblk)
  {
    ctx->cgetq = ctx->cgetv;
    last = 0;
    ioff = 0;
    nsblk = ctx->cgetq->cnblk;
  }
  else if (iblk == ctx->cgetv[1].cblk)
  {
    ctx->cgetq = &(ctx->cgetv)[1];
    last = 1;
    ioff = ctx->glob.maxblk;
    nsblk = ctx->cgetq->cnblk;
  }
  else if (last != 0)
  {
    ctx->cgetq = ctx->cgetv;
    nsblk = getqq(ctx, iblk, &(ctx->cgetq->cff), ctx->cgetq->cwt, ctx->ibkptr, ctx->ikmin, ctx->ivs);
    if (nsblk == 0)
    {
      *idgn = 0;
      return ctx->glob.nqn;
    }
    last = 0;
    ioff = 0;
    ctx->cgetq->cnblk = nsblk;
    ctx->cgetq->cblk = iblk;
    ctx->cgetq->csblk = 0;
  }
  else
  {
    ctx->cgetq = &(ctx->cgetv)[1];
    ioff = ctx->glob.maxblk;
    nsblk = getqq(ctx, iblk, &ctx->cgetq->cff, ctx->cgetq->cwt, &(ctx->ibkptr)[ioff],
                  &(ctx->ikmin)[ioff], &(ctx->ivs)[ioff]);
    if (nsblk == 0)
    {
      *idgn = 0;
      return ctx->glob.nqn;
    }
    last = 1;
    ctx->cgetq->cnblk = nsblk;
    ctx->cgetq->cblk = iblk;
    ctx->cgetq->csblk = ioff;
  }
  /*  check for request of size */
  if (indx <= 0)
  {
    *idgn = ctx->ibkptr[nsblk + ioff];
    return ctx->glob.nqn;
  }
  /*  search for sub-block */
  ix = indx - 1;
  ibgn = ctx->cgetq->csblk;
  iend = nsblk + ioff;
  if (ix < ctx->ibkptr[ibgn])
    ibgn = ioff;
  for (i = ibgn + 1; i < iend; ++i)
  {
    if (ix < ctx->ibkptr[i])
      break;
  }
  ctx->cgetq->csblk = ibgn = i - 1;
  /*  assemble quanta */
  ncod = ctx->ibkptr[i] - ctx->ibkptr[ibgn];
  ldgn = ctx->cgetq->cwt[0];
  ngsym = (short)setgsym(ctx, ctx->cgetq->cwt[3]);
  iqf = ctx->cgetq->cff;
  iqsp = getqs(ctx, ctx->ivs[ibgn], iqf, 0, 0, ctx->ixcom, ctx->iscom, &ivbase) - 1;
  iv = ctx->ixcom[XVIB];
  isym = ctx->ixcom[XSYM];
  n = ctx->ixcom[XNVAL];
  k = ix - ctx->ibkptr[ibgn];
  kq = ctx->ikmin[ibgn] + (k << 1);
  if (ngsym >= 3)
  {
    mgsym = (short)MOD(kq - ctx->ixcom[XLVAL] + ctx->ixcom[XISYM] + ngsym, ngsym);
    if (mgsym != 0)
    {
      if ((mgsym + mgsym) == ngsym)
      { /* B symmetry */
        ldgn = ctx->cgetq->cwt[2];
      }
      else
      { /* E symmetry */
        if (ctx->glob.esym)
          isym = 4;
        ldgn = ctx->cgetq->cwt[1];
      }
    }
  }
  lv = ctx->ixcom[XLVAL];
  if (lv != 0)
  { /* l-doubled state */
    if (ctx->glob.lsym)
      isym = 4;
    ncod = 1;
    if (kq == 0 && isym == 3)
    {
      if (iv == ivbase)
        ++iv;
      else
        --iv;
    }
  }
  if (ctx->glob.nqnn == 2)
  {
    if (ODD(isym))
      kq = -kq;
    ncod = 1;
  }
  else
  {
    if (isym <= 3)
      lower = (n + isym) & 1;
    else /* degenerate state */
      lower = (lv > 0) ? 1 : 0;
    if (ctx->glob.oblate)
    {
      iqn[2] = (short)kq;
      kq = n - kq + lower;
      ncod = -ncod;
    }
    else
    {
      iqn[2] = (short)(n - kq + lower);
    }
  }
  iqn[0] = (short)n;
  iqn[1] = (short)kq;
  k = ctx->glob.nqnn;
  if (ctx->glob.vibfmt)
    iqn[k++] = (short)iv;
  if (ctx->cgetq->cwt[4] > maxqn)
  {
    iqn[k++] = (short)iqsp;
    iqn[k] = (short)((iqf + 1) >> 1);
  }
  else
  {
    for (i = 0; k < maxqn; ++i)
    {
      iqn[k++] = (short)((ctx->iscom[i] + 1) >> 1);
      if (i == ctx->itsym && ngsym > 3)
        i += ngsym - 3;
    }
  }
  ldgn *= iqf + 1;
  if (ldgn < 0)
    ldgn = 0;
  *idgn = ldgn;
  return ncod;
} /* getqn */