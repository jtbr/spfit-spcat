/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 3 March 2004 */
/*   28 Sept.  2005: fix phase bug in getzitot for reduced matrix elements */
/* split from original spinv.c */

/* Functions dealing with spin algebra, spin coupling schemes, symmetry group handling, and statistical weights. */
/*
getsp(): Parses BCD spin codes and associates them with vibrational states.
setsp(): Finalizes the spin tables (spt arrays) with coupling information.
checksp(): Checks spin selection rules and calculates spin factors.
tensor(): Calculates spin coupling coefficients (6j, 9j symbols implicitly).
getll(): Determines tensor orders for spin couplings.
setgsym(): Sets global symmetry context (n-foldness, I<sub>tot</sub> parameters).
setwt(), getwt(): Sets/gets statistical weight for a given state symmetry and spin configuration.
testwt(): Tests if an operator can connect two states based on weights.
checkwt(): Compares two sets of weights.
blksym(): Calculates the D2 symmetry of an interaction block.
getmask(): Creates a selection mask for dircos based on K, L, and I<sub>tot</sub> symmetry.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calpgm.h"
#include "spinit.h"
#include "spinv_internal.h"
#include "SpinvContext.hpp"
#include "CalError.hpp"

/**
 * @brief Parses BCD encoded spin quantum numbers and associates them with a vibrational state.
 *        It manages a list of unique spin patterns (SSP structs) to save memory.
 * @param bcd_spin_codes Pointer to the BCD array containing spin codes. The first element [0]
 *                       holds flags and length. Subsequent elements encode 2*I values. Renamed ispnx.
 * @param pvinfo_vib_state Pointer to the SVIB structure for the vibrational state to associate
 *                         the spins with. Renamed pvinfo.
 * @return int Number of actual spin quantum numbers defined for this state (excluding N),
 *             adjusted for I_tot symmetry if applicable.
 */
int getsp(struct SpinvContext *ctx, const bcd_t *ispnx, SVIB *pvinfo)
{
  /*  subroutine to set up spin structure                        */
  /*     on entry:                                               */
  /*         ISPNX = spin code                                   */
  /*         IVIB  = vibrational quantum number                  */
  /*     on return:                                              */
  /*         JJS array of spin combinations                      */
  /*         (this array is initialized mostly in setsp)         */
  /*                                                             */
  SSP *ssp_now;
  short *jjs;
  size_t nl;
  int i, ii, jj, nspinv, idec, nspnx, nsstat, icmp, nset, isbig;
  int nitot, iret;
  short iis[MAXSPIN + 1];
  bcd_t itmp;

  jj = -1;
  nspnx = (int)(ispnx[0] & 0x7f);
  isbig = 0;
  itmp = ispnx[1];
  if ((int)(itmp & 0x0f) == 0)
  {
    jj = (int)(itmp >> 4);
    isbig = 1;
  }
  /*  decode spins */
  nspinv = 0;
  nl = 1;
  idec = isbig;
  for (i = 1; i <= MAXSPIN; ++i)
  {
    if (isbig == 0 && jj >= 0)
    {
      ii = jj;
      jj = -1;
    }
    else
    {
      ++idec;
      if (idec >= nspnx)
        break;
      itmp = ispnx[idec];
      ii = (int)(itmp & 0x0f);
      if (isbig != 0)
        ii = ii * 10 + jj;
      jj = (int)(itmp >> 4);
    }
    if (ii == 0)
      break;
    if (ii >= MAXII)
    {
      ii = 1;
    }
    else if (ii > 1)
    {
      nspinv = i;
    }
    nl *= (size_t)ii;
    --ii;
    iis[i] = (short)ii; /* TODO: should this be ii - 1? */
    if (ii > ctx->glob.mxspin)
      ctx->glob.mxspin = ii;
  }
  nsstat = (int)nl;
  i = (int)(nsstat << ctx->glob.msshft);
  if (i < 0 || (size_t)(i >> ctx->glob.msshft) != nl)
    throw InputError("spin problem too big", CalErrorCode::SpinDimensioning);
  nitot = 0;
  iret = nspinv;
  i = pvinfo->gsym;
  if (ODD(i))
  {
    nitot = i >> 1;
    if (nitot > nspinv)
    {
      iret += nitot;
      if (iret > MAXSPIN)
        iret = MAXSPIN;
      for (i = nspinv + 1; i <= iret; ++i)
        iis[i] = 0;
      nspinv = iret;
    }
    if (nitot > 3)
    {
      iret -= (nitot - 3);
    }
  }
  pvinfo->nspstat = nsstat;
  nset = nspinv + 1;
  iis[0] = (short)nset;
  ssp_now = &(ctx->ssp_head);
  do
  { /*  check previous combinations for match */
    icmp = 0;
    jjs = ssp_now->sspt;
    for (i = 0; i < nset; ++i)
    {
      icmp = iis[i] - jjs[i];
      if (icmp != 0)
        break;
    }
    if (icmp == 0)
      icmp = nitot - ssp_now->nitot;
    /*  match obtained, so return */
    if (icmp == 0)
    {
      pvinfo->spt = jjs;
      return iret;
    }
    ssp_now = ssp_now->next;
  } while (ssp_now != NULL);
  nl = (size_t)(nsstat + 1) * (size_t)nset * sizeof(short);
  jjs = (short *)calalloc(nl);
  nl = sizeof(SSP);
  ssp_now = (SSP *)calalloc(nl);
  ssp_now->next = ctx->ssp_head.next;
  ctx->ssp_head.next = ssp_now;
  ssp_now->sspt = jjs;
  ssp_now->ssize = nsstat;
  ssp_now->nitot = nitot;
  jjs[0] = (short)nset;
  for (i = 1; i < nset; ++i)
  {
    jjs[i] = iis[i];
  }
  pvinfo->spt = jjs;
  return iret;
} /* getsp */

/**
 * @brief Finalizes the spin tables (spt) for all unique spin patterns (SSP structs).
 *        For each unique spin pattern, it calculates and stores the coupling information:
 *        - Differences like 2*(N-F), 2*(J-F), 2*(F1-F), etc.
 *        - Minimum 2*F value for each spin combination.
 *        - I_tot symmetry components if applicable.
 *        This pre-computation is used by getqs to quickly retrieve quantum numbers.
 */
void setsp(struct SpinvContext *ctx)
{
  /*  subroutine to set up spin structure                        */
  /*   let: jjs = array of spin combinations in struct ssp.sspt  */
  /*         jjs[0]  = NSET = NSPIN + 1 = size of the set        */
  /*                   {N,J,F1..F} or {N,J..Itot,F}              */
  /*         jjs[1]..jjs[NSPIN] = S, I1, I2..                    */
  /*         jjs[k*NSET]..jjs[NSPIN-1+k*NSET] = 2*(N-F),2*(J-F)..*/
  /*         jjs[NSPIN+k*NSET] = minimum value of 2*F +(1)       */
  /*                                                             */
  SSP *ssp_now;
  ITMIX *pitmix0;
  EITMIX *pitmix;
  int nset, nspinv, ispsum, nn, iv, itot, ibgn, iend, ns, ns0, ifhalf;
  int minf0, minft, jj, ival, ibase, i, itmp, knt, nlim, nitot, ii;
  int itsym0, itptr0;
  short *iis, *jjs, *jjs0;
  ctx->nspin = 0;
  ssp_now = &ctx->ssp_head; /* head has no spins */
  while ((ssp_now = ssp_now->next) != NULL)
  { /* loop over spin sets */
    iis = ssp_now->sspt;
    nspinv = iis[0] - 1;
    for (iv = nspinv; iv > ctx->nspin; --iv)
    {
      if (iis[0] > 0)
      {
        ctx->nspin = iv;
        break;
      }
    }
  }
  ssp_now = &ctx->ssp_head; /* head has no spins */
  while ((ssp_now = ssp_now->next) != NULL)
  { /* loop over spin sets */
    iis = ssp_now->sspt;
    nset = iis[0];
    nspinv = nset - 1;
    if (nspinv > ctx->nspin)
      nspinv = ctx->nspin;
    ispsum = 0;
    for (iv = 1; iv <= nspinv; ++iv)
    {
      ispsum += iis[iv];
    }
    ifhalf = (ispsum & 1);
    iend = iis[nspinv];
    ival = iend + 1;
    ns = ival;
    pitmix0 = NULL;
    pitmix = NULL;
    nitot = ssp_now->nitot;
    if (nitot == 0)
    {
      nlim = nspinv;
      ibgn = iend;
      itsym0 = itptr0 = ctx->nspin + ctx->nspin + 1;
    }
    else
    { /* Itot coupling */
      itptr0 = ctx->nspin - 2;
      itsym0 = ctx->nspin - nitot;
      nlim = itsym0 + 1;
      ii = iend;
      for (i = nlim; i < nspinv; ++i)
      {
        if ((int)iis[i] != ii)
          break;
        ns *= ival;
      }
      if (i < ctx->nspin)
        throw InputError("spins under Itot should be the same", CalErrorCode::SpinDimensioning);
      pitmix0 = get_itmix(ii, nitot);
      if (nitot < 3)
        pitmix0 = NULL;
      iend = nitot * ii;
      ibgn = (iend & 1);
    }
    knt = ssp_now->ssize;
    ns0 = knt / ns;
    jjs = iis + nset;
    for (nn = -ispsum; nn <= ispsum; nn += 2)
    { /* select N quanta */
      minf0 = -ifhalf;
      if (nn < minf0)
        minf0 = nn; /* nn = 2 * (N - F) */
      if (pitmix0 != NULL)
        pitmix = pitmix0->mix;
      for (itot = ibgn; itot <= iend; itot += 2)
      {
        ns = ns0 * (itot + 1);
        if (pitmix != NULL && itot > ibgn)
          pitmix += 1;
        for (iv = 0; iv < ns; ++iv)
        { /* try all possible spins */
          jj = nn;
          jjs[0] = (short)nn;
          ival = iv;
          minft = minf0;
          for (i = 1; i < nlim; ++i)
          {
            ibase = iis[i] + 1;
            itmp = ival;
            ival /= ibase;
            itmp -= ival * ibase;
            jj += itmp - iis[i]; /* jj = J + N - S - 2F */
            if (jj < minft)
              minft = jj;
            jj += itmp; /* jj = 2 * (J - F) */
            jjs[i] = (short)jj;
          }
          jj += ival - itot;
          if (jj < minft)
            minft = jj;
          jj += ival;
          if (jj == 0)
          { /* this possibility gives F=0 */
            if (itptr0 < nspinv)
              jjs[itptr0 + 1] = (short)itot;
            jj = ifhalf - minft;
            minft = jj + (jj & 1);
            jjs[nspinv] = (short)minft;
            if (pitmix != NULL)
            {
              jjs[itsym0 + 1] = pitmix->qsym[0];
              for (i = itsym0 + 2; i <= itptr0; ++i)
                jjs[i] = 0;
              nitot = pitmix->n;
              jjs0 = jjs;
              for (i = 1; i < nitot; ++i)
              {
                jjs += nset;
                --knt;
                for (jj = 0; jj <= nspinv; ++jj)
                  jjs[jj] = jjs0[jj];
                jjs[itsym0 + 1] = pitmix->qsym[i];
              }
            }
            if (--knt <= 0)
            {
              ns0 = 0;
              break;
            }
            jjs += nset;
          }
        } /* end loop over all spins */
      } /* end loop over Itot */
    } /* end loop over N */
  } /* end loop over spin types */
} /* setsp */

/**
 * @brief Checks spin selection rules and calculates part of the spin-dependent factor.
 *        Also initializes I_tot related coefficients via setzitot for the first call for a spin set.
 * @param first_call_for_spin_set TRUE if this is the first time checksp is called for this specific
 *                                combination of spin indices (si1, si2) and vibrational states.
 *                                Used to initialize I_tot factors. Renamed first.
 * @param spin1_idx Index of the first spin involved in the operator.
 * @param spin2_idx Index of the second spin involved (0 if one-spin, >0 if two-spin).
 * @param vib_bra_spin_table Spin table (spt) for the bra vibrational state.
 * @param vib_ket_spin_table Spin table (spt) for the ket vibrational state.
 * @param matrix_element_factor_ptr Pointer to the matrix element factor, which will be multiplied by spin-related prefactors. Renamed zfac.
 * @return int 0 if spin selection rules are satisfied, otherwise a non-zero error code (2-6, or 1 for I_tot related issue).
 */
int checksp(struct SpinvContext *ctx, const BOOL first, int si1, int si2, const short *iiv1,
            const short *iiv2, double *zfac)
{
  int ii, iip;
  ii = iip = 0;
  if (si2 > 0)
  {
    if ((short)si2 < iiv1[0])
      ii = iiv1[si2];
    if ((short)si2 < iiv2[0])
      iip = iiv2[si2];
    if (ii == iip)
    {
      if (ii == 0)
        return 2;
      if (si2 == si1)
      {
        if (ii < 2)
          return 3;
        *zfac *= ctx->spfac2[ii];
      }
      *zfac *= ctx->spfac[ii];
    }
    else
    { /* ii != iip */
      if (ODD(ii + iip))
        return 4; /* check multiplicity */
      if (si2 > ctx->itsym)
        return 4;
    }
  }
  if (si2 != si1)
  {
    ii = iip = 0;
    if ((short)si1 < iiv1[0])
      ii = iiv1[si1];
    if ((short)si1 < iiv2[0])
      iip = iiv2[si1];
    if (ii == iip)
    {
      if (ii == 0)
        return 5;
      *zfac *= ctx->spfac[ii];
    }
    else
    { /* ii != iip */
      if (ODD(ii + iip))
        return 6; /* check multiplicity */
      if (si1 > ctx->itsym)
        return 6;
    }
  }
  if (first && si1 > ctx->itsym)
  {
    iip = 0;
    if (ctx->glob.nitot >= 3)
      iip = si1 - ctx->itsym;
    if (si1 == si2)
    {
      if (iip > 1)
        return 1;
      setzitot(ctx, 2, 0, 2, ii, ctx->glob.nitot); /* quadrupole */
    }
    else if (si2 > ctx->itsym)
    {
      if (iip > 2)
        return 1;
      setzitot(ctx, 1, 1, 0, ii, ctx->glob.nitot); /* 2-spin product */
      setzitot(ctx, 1, 1, 2, ii, ctx->glob.nitot);
    }
    else
    {
      if (iip > 1)
        return 1;
      setzitot(ctx, 1, 0, 1, ii, ctx->glob.nitot); /* 1-spin vector */
    }
  }
  return 0;
} /* checksp */

/**
 * @brief Calculates spin coupling factors (Wigner 6j, 9j symbols) for an operator.
 * @param matrix_element_factor_ptr Pointer to the matrix element factor to be modified by spin coeffs. Renamed zval.
 * @param iscom_bra_qns Array of 2*angular momentum values for bra state. Renamed iscom.
 * @param jscom_ket_qns Array of 2*angular momentum values for ket state. Renamed jscom.
 * @param lscom_tensor_orders Array of 2*tensor orders for intermediate and spin operators. Renamed lscom.
 * @param spin_index_map Maps operator spin index to actual spin index in iscom/jscom. Renamed smap.
 * @param num_coupled_pairs Number of spin coupling stages. Renamed npair.
 * @param alpha_sym_component I_tot alpha symmetry component. Renamed alpha.
 * @return int Always 0.
 */
int tensor(struct SpinvContext *ctx, double *zval, const int *iscom, const int *jscom,
           const int *lscom, const int *smap, int npair, int alpha)
{
  /*     subroutine to find spherical tensor spin coupling coefficients */
  /*     on entry: */
  /*         Z= initial multiplier */
  /*         NS= number of spins */
  /*         ISCOM,JSCOM= list of angular momenta *2 */
  /*         LLCOM= list of tensor orders *2 */
  /*         LSCOM = list of tensor orders for spins *2 */
  /*     on return: */
  /*         Z= modified multiplier */
  double zsq, z;
  int isgn, i, jjf, jji, nnf, llj, nni, lln, lls, ssf, ssi, ix, jj[9];

  if (npair <= 0)
    return 0;
  z = (*zval);
  isgn = 0;
  zsq = 1.;
  for (i = 0; i < npair; ++i)
  {
    ix = smap[i + i];
    lln = lscom[ix];
    nni = iscom[ix];
    nnf = jscom[ix];
    ix = smap[i + i + 1];
    lls = lscom[ix];
    ssi = iscom[ix];
    ssf = jscom[ix];
    if (alpha >= 0 && ix == ctx->itptr)
    {
      ix = ctx->itsym + ctx->nspin + 1;
      getzitot(&z, lls, iscom[ix], &lscom[ix],
               &iscom[ctx->itsym], &jscom[ctx->itsym], alpha, ctx->glob.nitot);
      i = ctx->nspin - 1;
    }
    llj = lscom[i];
    jji = iscom[i];
    jjf = jscom[i];
    if (lln != 0)
    {
      if (llj == 0)
      { /* T * U */
        isgn += jji + ssi + nnf;
        z *= c6jj(nni, lln, nnf, ssf, jji, ssi);
      }
      else if (lls == 0)
      { /* U = 1 */
        isgn += jjf + nni + ssi + llj;
        zsq *= jji + 1;
        zsq *= jjf + 1;
        z *= c6jj(nni, llj, nnf, jjf, ssi, jji);
      }
      else
      { /* T X U */
        isgn += 2;
        zsq *= 2.5;
        ix = (jjf + 1) * (llj + 1);
        zsq *= ix;
        zsq *= jji + 1;
        jj[0] = nnf;
        jj[1] = lln;
        jj[2] = nni;
        jj[3] = ssf;
        jj[4] = lls;
        jj[5] = ssi;
        jj[6] = jjf;
        jj[7] = llj;
        jj[8] = jji;
        z *= c9jj(jj);
      }
    }
    else if (lls != 0)
    { /* T = 1 */
      isgn += jji + nni + ssf + llj;
      zsq *= jji + 1;
      zsq *= jjf + 1;
      z *= c6jj(ssi, llj, ssf, jjf, nni, jji);
    }
  }
  if (ODD2(isgn))
    z = -z;
  *zval = z * sqrt(zsq);
  return 0;
} /* tensor */

/**
 * @brief Determines tensor orders for spin couplings and checks triangle inequalities.
 * @param llf_total_tensor_order 2 * Tensor order for the total F quantum number.
 * @param dir_cos_order Order of the direction cosine operator part.
 * @param n_tensor_order Order of the N-dependent tensor part.
 * @param k_delta_op Change in K for the operator.
 * @param spin1_idx Index of the first spin involved.
 * @param spin2_idx Index of the second spin involved.
 * @param lscom_tensor_orders_out Output array: lscom_out[0] to lscom_out[nspin-1] are 2*tensor orders for intermediate couplings (L_J, L_F1, ...).
 *                  lscom_out[nspin] to lscom_out[nspin+nspin] are 2*tensor orders for individual spin operators (L_S, L_I1, ...).
 * @param iscom_bra_qns Array of 2*angular momentum values for the bra state.
 * @param jscom_ket_qns Array of 2*angular momentum values for the ket state.
 * @return int Maximum spin index involved in the operator (up to nspin), or -1 if any triangle rule is violated.
 */
int getll(struct SpinvContext *ctx, const int llf_total_tensor_order, const int dir_cos_order,
          const int n_tensor_order, const int k_delta_op, const int spin1_idx, const int spin2_idx,
          int *lscom_tensor_orders_out, const int *iscom_bra_qns, const int *jscom_ket_qns)
{
  int *l_individual_spin_tensors_ptr;                                                                                                                                                                             /* Renamed lsscom */
  int i_spin_level, current_intermediate_L_times_2, ll_dir_cos_times_2, ll_n_tensor_times_2, max_L_N_operator_times_2, delta_2N_val, sum_2N_val, max_involved_spin_idx, error_return_flag, num_spins_minus_2_val; /* Renamed i, llj, lld, lln, llmax, jdif, jsum, maxspin, ierr, nm2 */

  error_return_flag = -1;
  max_L_N_operator_times_2 = 0;
  ll_dir_cos_times_2 = dir_cos_order << 1;
  ll_n_tensor_times_2 = n_tensor_order << 1;
  delta_2N_val = iscom_bra_qns[ctx->nspin] - jscom_ket_qns[ctx->nspin];
  if (delta_2N_val < 0)
    delta_2N_val = -delta_2N_val;

  if (ll_dir_cos_times_2 >= ll_n_tensor_times_2)
  {
    if (ll_n_tensor_times_2 < delta_2N_val)
      return error_return_flag;
    max_L_N_operator_times_2 = k_delta_op << 1;
    if (max_L_N_operator_times_2 < ll_dir_cos_times_2)
      max_L_N_operator_times_2 = ll_dir_cos_times_2;
  }
  else
  {
    if (ll_dir_cos_times_2 < delta_2N_val)
      return error_return_flag;
    max_L_N_operator_times_2 = k_delta_op << 1;
    if (max_L_N_operator_times_2 < ll_n_tensor_times_2)
      max_L_N_operator_times_2 = ll_n_tensor_times_2;
  }
  sum_2N_val = iscom_bra_qns[ctx->nspin] + jscom_ket_qns[ctx->nspin];
  if (max_L_N_operator_times_2 > sum_2N_val)
    return error_return_flag;

  if (ctx->nspin == 0)
  {
    lscom_tensor_orders_out[0] = ll_n_tensor_times_2;
    return 0;
  }

  max_involved_spin_idx = spin1_idx;
  if (spin1_idx > ctx->itsym && ctx->itsym < ctx->nspin)
    max_involved_spin_idx = ctx->nspin;

  l_individual_spin_tensors_ptr = &lscom_tensor_orders_out[ctx->nspin];
  l_individual_spin_tensors_ptr[0] = ll_n_tensor_times_2;
  for (i_spin_level = 1; i_spin_level <= ctx->nspin; ++i_spin_level)
  {
    l_individual_spin_tensors_ptr[i_spin_level] = 0;
  }
  if (spin1_idx > 0)
  {
    l_individual_spin_tensors_ptr[spin1_idx] = 2;
    if (spin2_idx > 0)
      l_individual_spin_tensors_ptr[spin2_idx] += 2;
  }

  current_intermediate_L_times_2 = l_individual_spin_tensors_ptr[0];
  num_spins_minus_2_val = ctx->nspin - 2;
  for (i_spin_level = 0; i_spin_level <= num_spins_minus_2_val; ++i_spin_level)
  {
    if (i_spin_level == ctx->itsym)
    {
      if (spin1_idx <= ctx->itsym)
      {
        current_intermediate_L_times_2 = 0;
      }
      else if (spin2_idx <= ctx->itsym || spin1_idx == spin2_idx)
      {
        current_intermediate_L_times_2 = l_individual_spin_tensors_ptr[spin1_idx];
      }
      i_spin_level = num_spins_minus_2_val;
    }
    else if (i_spin_level >= max_involved_spin_idx)
    {
      current_intermediate_L_times_2 = llf_total_tensor_order;
    }
    else
    {
      current_intermediate_L_times_2 -= l_individual_spin_tensors_ptr[i_spin_level + 1];
      if (current_intermediate_L_times_2 < 0)
        current_intermediate_L_times_2 = -current_intermediate_L_times_2;
    }
    lscom_tensor_orders_out[i_spin_level] = current_intermediate_L_times_2;
    delta_2N_val = iscom_bra_qns[i_spin_level] - jscom_ket_qns[i_spin_level];
    if (delta_2N_val < 0)
      delta_2N_val = -delta_2N_val;
    if (current_intermediate_L_times_2 < delta_2N_val)
      return error_return_flag;
    sum_2N_val = iscom_bra_qns[i_spin_level] + jscom_ket_qns[i_spin_level];
    if (current_intermediate_L_times_2 > sum_2N_val)
      return error_return_flag;
  }
  lscom_tensor_orders_out[ctx->nspin - 1] = llf_total_tensor_order;
  if (llf_total_tensor_order != 0)
    max_involved_spin_idx = ctx->nspin;

  return max_involved_spin_idx;
} /* getll */

/**
 * @brief Sets the global symmetry context (n-foldness, I_tot parameters) based on a gsym code.
 * @param group_symmetry_code The gsym code for a vibrational state (from SVIB struct).
 *                            Even values are standard symmetry, odd values indicate I_tot basis active.
 *                            Magnitude >> 1 gives n-foldness.
 * @return int The n-foldness of the symmetry group.
 */
int setgsym(struct SpinvContext *ctx, const int gsym)
{
  static int oldgsym = -1;
  static int nsym;
  int k, kk;
  if (gsym == oldgsym)
    return nsym;
  oldgsym = gsym;
  nsym = gsym >> 1;
  ctx->is_esym[0] = 0;
  for (k = 1, kk = nsym - 1; k < kk; ++k, --kk)
  {
    ctx->is_esym[k] = 1;
    ctx->is_esym[kk] = -1;
  }
  if (k == kk)
    ctx->is_esym[k] = 0;
  if (ODD(gsym))
  {
    ctx->itptr = ctx->nspin - 2;
    ctx->itsym = ctx->nspin - nsym;
    ctx->glob.nitot = nsym;
  }
  else
  {
    ctx->itptr = ctx->itsym = ctx->nspin + ctx->nspin + 1;
    ctx->glob.nitot = 0;
  }
  /* set up nominal spin coupling map */
  for (k = 0; k < ctx->nspin; ++k)
  {
    ctx->ismap[k + k] = k - 1; /* last N,J,F1 .. */
    if (k < ctx->itsym)
    {
      ctx->ismap[k + k + 1] = k + ctx->nspin + 1; /* S,I1,I2 .. */
    }
    else
    {
      ctx->ismap[k + k + 1] = ctx->itptr; /* Itot */
      break;
    }
  }
  ctx->ismap[0] = ctx->nspin; /* fix up position of N */
  return nsym;
} /* setgsym */

/**
 * @brief Sets statistical weights for a vibrational state based on symmetry options.
 * @param pvinfo_vib_state Pointer to the SVIB structure for the vibrational state.
 * @param vib_state_index_or_count If positive, it's the current vibrational state index.
 *                                 If negative, it's -N_total_vib_states, indicating to apply
 *                                 these weights to N_total_vib_states starting from pvinfo_vib_state.
 * @param axis_option_iax Internal axis code (1=z(a), 2=y(b), 3=x(c), 4=A_D2h_like, 5=B_D2h_like) for D2 or Cn symmetries.
 * @param weight_plus Statistical weight for "plus" or "even" parity states.
 * @param weight_minus Statistical weight for "minus" or "odd" parity states.
 * @param vibrational_symmetry_code_vsym Encodes vibrational symmetries for multiple states if vib_state_index_or_count is negative,
 *                                       or indicates if current state's weights should be swapped if positive.
 * @return int Always 0.
 */
int setwt(SVIB *pvinfov, const int ivib, const int iax,
          const int iwtpl, const int iwtmn, double vsym)
{
  /*  set weights */
  /*     on entry: */
  /*         IVIB= number of vibrational states * number of wang states */
  /*         IAX= axis of symmetry */
  /*         IWTPL= weight for even rotational states */
  /*         IWTMN= weight for odd rotational states */
  /*         VSYM= symmetry code for vibrational states */
  /*     on return: */
  /*         IWT= weights for vibration rotation states */
  SVIB *pvinfo;
  int n, nvsym, i, k;
  short ii, jj, itmp;

  n = 1;
  nvsym = 0;
  if (ivib < 0)
  {
    if (vsym > 0.5)
      nvsym = 1;
    n = -ivib;
  }
  pvinfo = pvinfov;
  for (i = 0; i < n; ++i)
  {
    ii = (short)iwtpl;
    jj = (short)iwtmn;
    if (nvsym != 0)
    {
      k = (int)(fmod(vsym, 10.) + 0.5);
      vsym = (vsym - k) / 10.;
      if (vsym < 0.)
        nvsym = 0;
      if (ODD(k))
      {
        itmp = ii;
        ii = jj;
        jj = itmp;
      }
    }
    switch (iax)
    {
    case 1: /*  WEIGHTS FOR Z */
      pvinfo->wt[0] = ii;
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = jj;
      pvinfo->wt[3] = ii;
      pvinfo->wt[4] = 2;
      ii = jj;
      break;
    case 2: /* WEIGHTS FOR B */
      pvinfo->wt[0] = ii;
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = ii;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 3;
      break;
    case 3: /* WEIGHTS FOR X */
      pvinfo->wt[0] = ii;
      pvinfo->wt[1] = ii;
      pvinfo->wt[2] = jj;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 2;
      break;
    case 4: /* A1, A2 only */
      pvinfo->wt[0] = ii;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 3;
      break;
    case 5: /* B1, B2 only */
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = ii;
      pvinfo->wt[4] = 3;
      break;
    } /* end switch */
    ++pvinfo;
  }
  return 0;
} /* setwt */

/**
 * @brief Gets the appropriate statistical weight for a given state symmetry and spin configuration.
 * @param pvinfo_vib_state Pointer to the SVIB structure for the vibrational state.
 * @param d2_symmetry_idx D2 symmetry index (0-3) of the rotational state.
 * @param spin_state_idx If 0, requests general check for existence of states for this d2_symmetry_idx.
 *                       If >0, it's the spin pattern index (1-based) from the spin table (spt).
 *                       If <0, it's -alpha_symmetry_code for I_tot E-state weight lookup.
 * @param weights_out Output array of 3 integers:
 *                    weights_out[0] = weight for A/B type rotational states (or primary component).
 *                    weights_out[1] = weight for E type rotational states (first E component, e.g. E_a).
 *                    weights_out[2] = weight for other E or B type states (e.g. E_b or B for D4/D6).
 *                    Unused components are set to -1.
 * @return int The primary statistical weight (weights_out[0]), or sum of weights if I_tot basis and multiple components allowed.
 *             Returns <=0 if no states/weight for this symmetry/spin.
 */
int getwt(struct SpinvContext *ctx, SVIB *pvinfo, const int isym, const int iispin, int *ivwt)
{
  int iwt, jsym, nset, k, msym, nsym;
  short *jjs;
  if (iispin == 0)
  {
    iwt = 0;
    if (pvinfo->knmax >= pvinfo->knmin[isym] ||
        pvinfo->knmax >= pvinfo->knmin[3 - isym])
      iwt = 1;
    return iwt;
  }
  k = 0;
  jsym = isym & 3;
  iwt = pvinfo->wt[4];
  if (iwt != 0)
  {
    /* resolve Itot symmetries */
    if (iispin > 0)
    {
      jjs = pvinfo->spt;
      nset = jjs[0];
      k = jjs[iispin * nset + ctx->itsym + 1];
    }
    else
    {
      k = -iispin;
    }
    if (iwt > 0 && ODD2(k))
      jsym ^= iwt;
  }
  ivwt[0] = pvinfo->wt[jsym];
  ivwt[1] = ivwt[2] = -1;
  nsym = pvinfo->gsym >> 1;
  if (nsym < 3)
    return ivwt[0];
  if (ODD(nsym))
  {                           /* nsym = 3, 5 */
    ivwt[1] = pvinfo->ewt[0]; /* E state */
  }
  else
  { /* nsym = 4, 6 */
    msym = ctx->isoddk[jsym] - (int)(pvinfo->lvqn);
    if (iwt != 0)
      msym += k >> 2;
    if (ODD(msym))
    {
      ivwt[0] = -1;             /* no A state */
      ivwt[1] = pvinfo->ewt[1]; /* E state */
      if (nsym == 6)
        ivwt[2] = pvinfo->wt[jsym]; /* B state */
    }
    else
    {
      ivwt[1] = pvinfo->ewt[0]; /* E state */
      if (nsym == 4)
        ivwt[2] = pvinfo->wt[jsym ^ 2]; /* B state */
    }
  }
  iwt = ivwt[0];
  if (pvinfo->gsym > 6)
  {
    if (iwt <= 0)
      iwt = ivwt[1];
    if (iwt <= 0)
      iwt = ivwt[2];
  }
  return iwt;
} /* getwt */

/**
 * @brief Tests if an operator can connect two vibrational states based on their statistical weights
 *        and the operator's symmetry.
 * @param pvib1_bra Pointer to SVIB structure for the bra state.
 * @param pvib2_ket Pointer to SVIB structure for the ket state.
 * @param operator_d2_symmetry D2 symmetry of the operator (0-3).
 * @param operator_alpha_symmetry I_tot alpha symmetry component of the operator.
 * @return BOOL TRUE if the connection is forbidden by weight/symmetry rules, FALSE otherwise.
 */
BOOL testwt(struct SpinvContext *ctx, SVIB *pvib1, SVIB *pvib2, int isym, int alpha)
{ /* return true if all weight sets do not match */
  int ii, jj, kk, ix, jx, iwt4, iwt[3], jwt[3];
  if (pvib1->gsym != pvib2->gsym)
    return TRUE;
  iwt4 = pvib1->wt[4];
  if ((short)iwt4 != pvib2->wt[4])
    return TRUE;
  ix = -1;
  jx = ix - (alpha << 2);
  for (kk = 0; kk < 2; ++kk)
  {
    for (ii = 0; ii < 4; ++ii)
    {
      getwt(ctx, pvib1, ii, ix, iwt);
      jj = ii ^ isym;
      getwt(ctx, pvib2, jj, jx, jwt);
      if (checkwt(ctx, iwt, jwt) == 0)
        return FALSE;
    }
    if (iwt4 <= 0)
      break;
    jx -= 2;
  }
  return TRUE;
} /* testwt */

/**
 * @brief Checks if any component of two sets of statistical weights are equal and positive.
 *        Used by testwt to determine if an operator can connect states.
 * @param weights_bra Array of 3 weights for the bra state (A/B, E1, E2/B').
 * @param weights_ket Array of 3 weights for the ket state.
 * @return int 0 if any pair of corresponding positive weights are equal, 1 otherwise.
 */
int checkwt(struct SpinvContext *ctx, int *iwt, int *jwt)
{
  /* return 0 if any pair of weights are equal */
  int ii, jj, nn;
  nn = ctx->glob.maxwt;
  if (nn == 0)
    return (iwt[0] - jwt[0]);
  for (ii = nn; ii >= 0; --ii)
  {
    if (iwt[ii] <= 0)
      continue;
    for (jj = nn; jj >= 0; --jj)
    {
      if (iwt[ii] == jwt[jj])
        return 0;
    }
  }
  return 1;
} /* checkwt */

/**
 * @brief Calculates the overall symmetry of an interaction block.
 * @param ixcom Quantum number block for the bra state.
 * @param jxcom Quantum number block for the ket state.
 * @return unsigned int The symmetry of the interaction (0=A, 1=Bx, 2=By, 3=Bz in D2),
 *                      result of XORing bra and ket state symmetries.
 */
unsigned int blksym(const int *ixcom, const int *jxcom)
{ /* get block symmetry from state symmetry */            /* Original comment */
  return (unsigned int)(3 & (ixcom[XSYM] ^ jxcom[XSYM])); /* XOR symmetries and mask to 2 bits */
} /* blksym */

/**
 * @brief Creates a mask for dircos to select allowed K and L transitions based on symmetry and operator type.
 * @param xbra_qns Quantum number block for the bra state. Renamed xbra.
 * @param xket_qns Quantum number block for the ket state. Renamed xket.
 * @param k_delta_op Total Delta K for the operator. Renamed kd.
 * @param l_delta_op Change in l quantum number (Delta l) for the operator. Renamed ldel.
 * @param operator_flags Operator property flags (MODD, MSYM2, MLZ, MIDEN from SPAR/DIRCOS). Renamed loff.
 * @param alpha_sym_component I_tot alpha symmetry component (0 if not I_tot or alpha=0 component). Renamed alpha.
 * @return int Mask value: bit 0 for (K'=K+DK, L'=L+DL), bit 1 for (K'=DK-K, L'=DL-L), bit 2 for (K'=K-DK, L'=L-DL).
 *         Returns 0 if no transitions are allowed.
 */
int getmask(struct SpinvContext *ctx, const int *xbra, const int *xket,
            const int kd, const int ldel,
            const int loff, const int alpha)
{
  /* CREATE MASK FOR DIRCOS */
  /* KD IS TOTAL K QUANTIUM CHANGE */
  int ldif, lbra, lket, mask, kbra, kket, nitot, kk, ksbra, ksket;
  mask = 7;
  nitot = ctx->glob.nitot;
  if (nitot >= 3)
  { /* check Itot symmetry */
    ksbra = xbra[XISYM];
    ksket = xket[XISYM];
    kk = (xbra[XIQN] ^ xket[XIQN]) >> 1;
    if (TEST(loff, MIDEN))
    { /* Identity operator */
      if (kk != 0)
        return 0;
    }
    else if (EVEN(kk))
    { /* A1|A1, A1|E, E|E, A2|A2, ..*/
      if (TEST(loff, MSYM2))
        return 0;
    }
    else
    { /* kk is ODD */
      if (TEST(loff, MSYM2))
        return 0; /* <A2| I_alpha |E> */
      if (ctx->is_esym[ksbra] == 0 && ctx->is_esym[ksket] == 0)
        return 0; /*  <A1 |I_alpha| A2> == 0 */
    }
    mask = 0;
    kk = ksket - ksbra;
    if (alpha == 0)
    {
      if (kk == 0)
        mask = 5;
    }
    else if (alpha < nitot)
    {
      if (MOD(kk + alpha, nitot) == 0)
        mask = 1;
      if (MOD(kk - alpha, nitot) == 0)
        mask += 4;
    }
    else
    { /* quasi-diagonal */
      if (MOD(ksbra + ksket - alpha, nitot) == 0)
        mask = 2;
    }
    if (mask == 0)
      return mask;
  }
  /* mask bit 0 (val = 1): kbra = kket + kd, lbra = lket + ldel */
  /* mask bit 2 (val = 4): kbra = kket - kd, lbra = lket - ldel */
  /* mask bit 1 (val = 2): kbra = kd - kket, lbra = ldel - lket */
  lbra = xbra[XLVAL];
  lket = xket[XLVAL];
  ldif = lbra - lket;
  if (ldel == 0)
  { /* operator diagonal in l */
    if (ldif != 0)
      mask &= 2; /* clear bits 0 and 2 */
  }
  else
  { /* operator off-diagonal in l */
    if (ldif == 0)
    {
      mask &= 2; /* clear bits 0 and 2 */
    }
    else if (ldif != ldel)
    {
      mask &= 6; /* clear bit 0 */
    }
    else if (-ldif != ldel)
    {
      mask &= 3; /* clear bit 2 */
    }
  }
  if ((mask & 2) == 0)
    return mask; /* bit 1 not set */
  mask &= 5;     /* clear bit 1 temporarily */
  /* quasi-diagonal, KD = KBRA + KKET */
  if (lbra + lket != ldel)
    return mask;
  kbra = xbra[XKBGN];
  kket = xket[XKBGN];
  kk = kd - kbra - kket;
  if (kk < 0)
    return mask;
  if (EVEN(xket[XIQN]))
    kk -= 2;
  if (EVEN(xbra[XIQN]))
    kk -= 2;
  if (kk < 0)
    return mask;
  return (mask + 2);
} /* getmask */
