/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 3 March 2004 */
/*   28 Sept.  2005: fix phase bug in getzitot for reduced matrix elements */
/* split from original spinv.c */

/* Functions for linear algebra helpers, sorting eigenvalues/eigenvectors, and specific Hamiltonian ordering schemes. */
/*
ordham(): Orders eigenvalues/H-diagonals within sub-blocks.
fixham(): Applies a permutation to reorder eigenvalues/vectors.
kroll(): Checks for and corrects Hamiltonian K-roll-over.
bestk(): Performs K-based sorting.
dclr(): Clears an array (very generic, but often used with matrices).
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

/**
 * @brief Orders eigenvalues within sub-blocks to roughly match the diagonal of the Hamiltonian.
 *        Used when glob.idiag = 2 or 5 for initial ordering before fixham.
 * @param nn Dimension of the Hamiltonian.
 * @param mask Eigenvector mixing coefficient or projection quantum number (not directly used for sorting here, but passed).
 * @param egy Array of eigenvalues (sorted by energy within sub-blocks initially by hdiag),
 *            also contains original diagonal elements of H if idiag=2,5. This function sorts 'egy' based on its values.
 * @param isblk Array of pointers to the start of each sub-block (Wang or vibrational).
 * @param iswap Output array storing the permutation required to achieve the desired order.
 * @return int Always 0.
 */
int ordham(const int nn, short *mask, double *egy,
           const short *isblk, short *iswap)
{
  /* subroutine to order eigenvalues within sub-block like diagonal of */                   /* Original comment */
  /*      Hamiltonian */                                                                    /* Original comment */
  double temp_egy_val;                                                                      /* Renamed tmp */
  int current_sub_block_idx, i_row, i_check, i_min_idx, sub_block_end_idx, inext_isblk_val; /* Renamed iblk, i, ii, iq, is, inext */
  short mcmp_val;                                                                           /* Renamed mcmp */

  current_sub_block_idx = 1;  /* sub-block counter, seems 1-based from isblk[1] */
  i_row = 0;                  /* Current row/state index we are trying to place */
  inext_isblk_val = isblk[1]; /* Index of the start of the *next* sub-block */
  for (sub_block_end_idx = 1; sub_block_end_idx < nn; ++sub_block_end_idx)
  { /* Loop through states, up to nn-1 */
    if (sub_block_end_idx == inext_isblk_val)
    {                                                   /* If we reached the end of the current sub-block */
      inext_isblk_val = isblk[++current_sub_block_idx]; /* Get end of next sub-block */
      iswap[i_row] = (short)i_row;                      /* The last element of a sub-block is already in place relative to itself */
    }
    else
    {                            /* Within a sub-block, find element with minimum 'egy' value */
      temp_egy_val = egy[i_row]; /* Current minimum 'egy' value */
      i_min_idx = i_row;         /* Index of that minimum */
      mcmp_val = mask[i_row];    /* Get mask (iqnsep) value for current state */
      for (i_check = sub_block_end_idx; i_check < inext_isblk_val; ++i_check)
      { /* Search remaining elements in this sub-block */
        if (egy[i_check] < temp_egy_val && mask[i_check] == mcmp_val)
        { /* If smaller 'egy' found and same mask type */
          temp_egy_val = egy[i_check];
          i_min_idx = i_check;
        }
      }
      iswap[i_row] = (short)i_min_idx; /* Record where the element for row 'i_row' was found */
      if (i_min_idx > i_row)
      { /* If the minimum was not at 'i_row', swap it into 'i_row' */
        egy[i_min_idx] = egy[i_row];
        egy[i_row] = temp_egy_val;
      }
    }
    i_row = sub_block_end_idx; /* Move to next row to process */
  }
  /* The last element (nn-1) needs its iswap value set */
  if (nn > 0)
    iswap[nn - 1] = (short)(nn - 1); /* Assuming it's correctly placed by previous logic or only one left */

  return 0;
} /* ordham */

/**
 * @brief Reorders eigenvalues, eigenvectors, and mixing coefficients based on a permutation vector.
 * @param ndm Leading dimension of matrix t (eigenvectors).
 * @param nn Dimension of the Hamiltonian / number of states.
 * @param t Eigenvector matrix (modified in place).
 * @param egy Eigenvalue array (modified in place).
 * @param p Mixing coefficient array (pmix) (modified in place).
 * @param iswap Permutation array from ordham or bestk. iswap[i] is the original index of the state that should go to position i.
 * @return int Always 0.
 */
int fixham(const int ndm, const int nn, double *t, double *egy, double *p, const short *iswap)
{
  /* subroutine to order eigenvalues within sub-block like diagonal of */ /* Original comment (slightly misleading, it applies a pre-calculated order) */
  /*      Hamiltonian using permutation found with ordham */              /* Original comment */
  int i_current, original_idx_for_current_pos;                            /* Renamed i, iq */
  int loop_end_val;                                                       /* Renamed is */

  loop_end_val = nn - 2; /* Loop from nn-2 down to 0 */
  if (loop_end_val >= 0)
  {
    for (i_current = loop_end_val; i_current >= 0; --i_current)
    {
      original_idx_for_current_pos = iswap[i_current]; /* Get original index of the element that should be at i_current */
      if (original_idx_for_current_pos > i_current)
      { /* If it's not already in place, swap it */
        /* etswap swaps eigenvectors, eigenvalues, and pmix values between original_idx_for_current_pos and i_current */
        etswap(ndm, nn, original_idx_for_current_pos, i_current, t, egy, p);
      }
    }
  }
  return 0;
} /* fixham */

/**
 * @brief Checks for Hamiltonian roll-over (diagonal elements not monotonically increasing/decreasing).
 *        If roll-over is detected, it modifies the Hamiltonian to enforce monotonicity
 *        and zeros out corresponding off-diagonal elements.
 * @param hamiltonian_dim Leading dimension of matrix t (Hamiltonian). Renamed nsizd.
 * @param t Hamiltonian matrix (modified in place if roll-over).
 * @param num_sub_blocks Number of sub-blocks. Renamed nsblk.
 * @param sub_block_pointers Array of pointers to the start of each sub-block. Renamed sbkptr.
 * @param kmin_values Array of minimum K values for each sub-block. Renamed kmin.
 * @return BOOL TRUE if roll-over was detected and corrected, FALSE otherwise.
 */
BOOL kroll(const int nsizd, double *t, const int nsblk, const short *sbkptr, const short *kmin)
{
  /* subroutine to make sure that diagonal elements of the Hamiltonian */
  /*    are monotonically increasing (decreasing for oblate) */
  /*    after K=KTHRSH */
  double vall, tlast, tmp, *ptmp;
  long ndmt;
  int ibgn, iend, i, k, n, ixx;
  BOOL roll;
  const double zero = 0.0;

  ndmt = nsizd + 1;
  roll = FALSE;
  for (ixx = 0; ixx < nsblk; ++ixx)
  {
    i = 5 - kmin[ixx];
    ibgn = sbkptr[ixx];
    if (i >= 2)
      ibgn += i >> 1;
    iend = sbkptr[ixx + 1] - 1;
    if (ibgn >= iend)
      continue;
    ptmp = &t[ibgn * ndmt];
    tmp = tlast = *ptmp;
    ++ibgn;
    for (i = ibgn; i <= iend; ++i)
    {
      ptmp += ndmt;
      tlast = tmp;
      tmp = (*ptmp);
      if (tmp < tlast)
        break;
    }
    if (i > iend)
      continue;
    k = kmin[ixx] + ((i - sbkptr[ixx]) << 1);
    printf(" roll-over at K = %d\n", k);
    roll = TRUE;
    vall = 1e15;
    for (; i <= iend; ++i)
    {
      n = i - 1;
      dcopy(n, &zero, 0, &t[i], nsizd);
      n = nsizd - i;
      memset(&t[i * ndmt], 0, n * sizeof(double));
      tlast += vall;
      *ptmp = tlast;
      ptmp += ndmt;
    }
  }
  return roll;
} /* kroll */

/**
 * @brief Performs K-based sorting of eigenvalues and eigenvectors, calculating K expectation values.
 * @param ndm Leading dimension of matrix t.
 * @param nsize Dimension of the Hamiltonian.
 * @param iqnsep Eigenvector mixing coefficients or projection labels from hdiag.
 * @param ibkptr Array of pointers to the start of each vibrational sub-block.
 * @param itau Array of tau-like quantum numbers (Ka-Kc proxy) or K values for each basis state.
 * @param iswap Output array for the permutation vector.
 * @param t Eigenvector matrix (modified in place).
 * @param egy Eigenvalue array (modified in place).
 * @param pmix Mixing coefficient array (modified in place, will store <K_op^2> / <sum |c_i|^2>).
 * @param wk Work array.
 * @return int Always 0.
 */
int bestk(const int ndm, const int nsize, short *iqnsep, short *ibkptr,
          short *itau, short *iswap, double *t, double *egy, double *pmix, double *wk)
{
#define MAXNK 4
  double ele;
  long ndml, tindx, tindx0;
  int i, iz, kd, k, ii, jj, ibgn, iend, iblk, is, iq, nk;
  short mcmp;
  ndml = ndm;
  memset(wk, 0, sizeof(double) * nsize);
  for (i = 0; (ibgn = ibkptr[i]) < nsize; ++i)
  {
    iend = ibkptr[i + 1] - 1;
    tindx0 = ibgn * (ndml + 1);
    for (ii = ibgn; ii <= iend; ++ii)
    {
      /* calculate expectation for PzPz */
      kd = (int)(itau[ii] >> 2);
      tindx = tindx0;
      for (jj = ibgn; jj <= iend; ++jj)
      {
        ele = kd * t[tindx];
        wk[jj] += ele * ele;
        tindx += ndml;
      }
      ++tindx0;
    }
  }
  for (i = 0; i < nsize; ++i)
  {
    wk[i] = sqrt(wk[i] / pmix[i]);
  }
  iend = ibkptr[1];
  iblk = 1;
  i = 0;
  for (is = 1; is < nsize; ++is)
  {
    if (is == iend)
    {
      iend = ibkptr[++iblk];
    }
    else
    { /* find min wk and itau */
      mcmp = iqnsep[i];
      iq = i;
      ele = wk[i];
      iz = i;
      jj = itau[i];
      for (ii = is; ii < iend; ++ii)
      {
        if (mcmp != iqnsep[ii])
          continue;
        if (wk[ii] < ele)
        {
          ele = wk[ii];
          iq = ii;
        }
        if (itau[ii] < (short)jj)
        {
          iz = ii;
          jj = itau[ii];
        }
      }
      if (iq > i)
      {
        wk[iq] = wk[i];
        wk[iq] = ele;
        etswap(ndm, nsize, iq, i, t, egy, pmix);
      }
      iswap[i] = (short)iz;
      if (iz > i)
      {
        itau[iz] = itau[i];
        itau[i] = (short)jj;
      }
    }
    i = is;
  }
  iend = ibkptr[1];
  iblk = 1;
  i = 0;
  nk = 0;
  for (is = 1; is < nsize; ++is)
  {
    if (is == iend)
    {
      iend = ibkptr[++iblk];
      nk = 0;
    }
    else
    {
      if (nk <= 0)
      {
        kd = itau[i] >> 2;
        nk = 1;
        for (ii = is; ii < iend; ++ii)
        {
          k = itau[ii] >> 2;
          if (k > kd)
            break;
          /* find basis with K == kd, nk = the number of basis with same K */
          ++nk;
        } /* loop ii */
      }
      if (nk > 1)
      {
        mcmp = iqnsep[i];
        ele = egy[i];
        iq = i;
        for (iz = 1; iz < nk; ++iz)
        {
          ii = iz + i;
          if (egy[ii] < ele && mcmp == iqnsep[ii])
          {
            ele = egy[ii];
            iq = ii;
          }
        }
        if (iq > i)
        {
          ele = wk[iq];
          wk[iq] = wk[i];
          wk[iq] = ele;
          etswap(ndm, nsize, iq, i, t, egy, pmix);
        }
      }
      --nk;
    }
    i = is;
  }
  return 0;
} /* bestk */

/**
 * @brief Clears a 2D array (or a block of it) by setting elements to zero.
 * @param n1 Number of rows (or elements if ix=0 or vec is 1D conceptually).
 * @param n2 Number of columns (if ix=1, effectively n1*n2 total elements).
 * @param vec Pointer to the array/vector to be cleared.
 * @param ix Stride/increment: 0 for contiguous, 1 for column-major like access if n1 is leading dim.
 *           (dcopy behavior: if ix=0, all n1 elements are zeroed. if ix=1, elements vec[0], vec[1]... zeroed.
 *            If this is meant for a matrix (N1 rows, N2 cols) with ix=1 for column-major,
 *            it would clear the first N1*N2 elements contiguously.
 *            However, dcopy(count, src, inc_src, dest, inc_dest) means this is dcopy(N, &zero, 0, pvec, ix).
 *            If ix=1, it clears N elements contiguously. If ix=0, it's not standard for dcopy.
 *            Given dcopy(&zero,0,...), src value is always zero, inc_src is irrelevant.
 *            This function clears n1*n2 elements of vec, advancing by ix for each of n1*n2 elements.
 *            This seems to be intended to clear an n1*n2 matrix when ix=1 means contiguous.
 *            If ix > 1, it would clear strided.
 *            The usage in hamx is dclr(nsize, nsize, t, 1) and dclr(nsize, n, dedp, 1),
 *            where 1 means contiguous clearing of nsize*nsize or nsize*n elements.
 *            This dclr itself calls dcopy(count, &zero, 0, pvec, ix).
 *            If ix=1 in dclr, then dcopy uses incdest=1. If ix=0 in dclr, dcopy uses incdest=0 (setting only pvec[0]).
 *            The hamx calls use ix=1, so it's contiguous clearing.
 *            The 'ix' parameter in dclr seems to be directly passed as incDest to dcopy.
 *            So dclr(N, M, vec, 1) clears N*M elements of vec contiguously.
 *            dclr(N, M, vec, stride) clears N*M elements with specified stride.
 *            Actually, dcopy takes (count, src_ptr, src_inc, dest_ptr, dest_inc).
 *            Here, dcopy((int)nsq, &zero, 0, pvec, ix) will set (int)nsq elements of pvec (starting at pvec[0])
 *            to zero, where each element is pvec[k*ix]. This seems to be clearing with a stride 'ix'.
 *            In hamx, ix=1 is used, meaning contiguous clearing.
 * )
 * @return int Always 0.
 */
int dclr(const int n1, const int n2, double *vec, const int ix)
{ /*  clear a N1*N2 block */            /* Original comment */
  static long nbig_chunk_size = 0x7ff0; /* Max elements for dcopy at once (approx 32k, for 16-bit int limit?) */
  long n_elements_to_clear;             /* Renamed nsq */
  double *current_vec_ptr;              /* Renamed pvec */
  const double zero = 0.0;

  current_vec_ptr = vec;
  n_elements_to_clear = n1 * (long)n2;          /* Total number of elements to clear */
  while (n_elements_to_clear > nbig_chunk_size) /* If total exceeds chunk size for dcopy */
  {
    dcopy((int)nbig_chunk_size, &zero, 0, current_vec_ptr, ix);
    current_vec_ptr += nbig_chunk_size; /* <--- Ensure this simple advance is used */
    n_elements_to_clear -= nbig_chunk_size;
  }
  dcopy((int)n_elements_to_clear, &zero, 0, current_vec_ptr, ix); /* Clear remaining elements */
  return 0;
} /* dclr */