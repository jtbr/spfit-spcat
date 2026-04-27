/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CalCat.hpp"
#include "splib/calpgm.h"
#include "common/CalError.hpp"

int CalCat::qnfmt(short *iqu, int nqn, char *sqn)
{ /* subroutine to do quick conversion of quantum integers to characters */
  static int czero = (int)'0';
  int i, k, ix;
  char tqn1, tqn2;

  for (k = 0; k < nqn; ++k) {
    i = iqu[k];
    if (i < 0) {
      i = -i;
      ix = i / 10;
      i -= ix * 10;
      tqn1 = (char)(i + czero);
      if (ix == 0) {
        tqn2 = '-';
      } else if (ix < 27) {
        tqn2 = (char)(ix + ((int)'a' - 1));
      } else {
        tqn1 = tqn2 = '*';
      }
    } else {
      ix = i / 10;
      i -= ix * 10;
      tqn1 = (char)(i + czero);
      if (ix == 0) {
        tqn2 = ' ';
      } else if (ix < 10) {
        tqn2 = (char)(ix + czero);
      } else if (ix < 36) {
        tqn2 = (char)(ix + ((int)'A' - 10));
      } else {
        tqn1 = tqn2 = '*';
      }
    }
    ix = k + k;
    sqn[ix] = tqn2;
    sqn[ix + 1] = tqn1;
  }
  return 0;
}

int CalCat::simt(int isiz, int jsiz, double *s, double *t, double *u, double *wk)
{
  int i, j;
  double *sij, *si, *tc, *uc;
  /*    subroutine to do similarity transform */
  /*      S' = transpose(T) * S * U */
  if (isiz > 1) { /* do transform on the left */
    sij = s;
    for (j = 0; j < jsiz; ++j) {
      dcopy(isiz, sij, 1, wk, 1);
      tc = t;
      for (i = 0; i < isiz; ++i) {
        *sij = ddot(isiz, tc, 1, wk, 1);
        tc += isiz;
        ++sij;
      }
    }
  }
  if (jsiz > 1 && u != NULL) { /* do transform on the right */
    si = s;
    for (i = 0; i < isiz; ++i) {
      dcopy(jsiz, si, isiz, wk, 1);
      sij = si;
      uc = u;
      for (j = 0; j < jsiz; ++j) {
        *sij = ddot(jsiz, wk, 1, uc, 1);
        uc += jsiz;
        sij += isiz;
      }
      ++si;
    }
  }
  return 0;
}

SBLK *CalCat::ibufof(const int ipos, const unsigned int ndel, SBLK *blk)
{
  SBLK *pmblk, *pdblk;
  double *pvec;
  long ldisk, len, n, lret;
  int k;

  if (ipos >= 0) {
    pmblk = &blk[ipos];
    if (pmblk->egyblk != NULL)
      return pmblk;
  }
  k = m_cache.mempos;
  m_cache.mempos = 0;
  pmblk = blk;
  while (pmblk->ixblk > 0) { /* find first empty block */
    ++m_cache.mempos;
    ++pmblk;
  }
  if (pmblk->egyblk == NULL) { /* write memory contents to disk */
    if (m_cache.maxrec == 0) {
      m_cache.scratch = tmpfile();
      if (m_cache.scratch == NULL)
        throw NumericError("scratch file open error", CalErrorCode::ScratchFileIo);
      m_cache.nbsav = m_cache.mempos;
      m_cache.maxdm = pmblk->nsizblk;
      ldisk = (long)ndel;
      if (blk->eigblk != NULL)
        ldisk += (long)m_cache.maxdm;
      m_cache.lsizb = (long)m_cache.maxdm * ldisk * (int)sizeof(double);
    }
    ldisk = (long)(m_cache.mempos - m_cache.nbsav);
    if (ldisk != 0)
      ldisk *= m_cache.lsizb;
    fseek(m_cache.scratch, ldisk, SEEK_SET);
    m_cache.mempos = k + 1;
    if (m_cache.mempos == m_cache.orgpos)
      ++m_cache.mempos;
    if (m_cache.mempos >= m_cache.nbsav)
      m_cache.mempos = (m_cache.orgpos != 0) ? 0 : 1;
    pdblk = pmblk;
    pmblk = &blk[m_cache.mempos];
    len = (long)ndel;
    n = (long)pmblk->nsizblk;
    len *= n;
    pvec = pmblk->egyblk;
    lret = 0;
    if (pvec != NULL) {
      lret = (long)fwrite(pvec, sizeof(double), (size_t)len, m_cache.scratch);
      if (lret == len && pmblk->eigblk != NULL) {
        if (ldisk == m_cache.maxrec)
          n = (long)m_cache.maxdm;
        len = n * n;
        lret = (long)fwrite(pmblk->eigblk, sizeof(double), (size_t)len,
                            m_cache.scratch);
      }
      if (ldisk == m_cache.maxrec && lret == len) {
        /* fill in gaps with junk FIRST time */
        n = (long)m_cache.maxdm - (long)pmblk->nsizblk;
        len = (long)ndel;
        if (n > 0) {
          len *= n;
          lret = (long)fwrite(pvec, sizeof(double), (size_t)len, m_cache.scratch);
        }
        m_cache.maxrec += m_cache.lsizb;
      }
    }
    if (lret != len)
      throw NumericError("heap write error", CalErrorCode::ScratchFileIo);
    pdblk->nsizblk = pmblk->nsizblk;
    pdblk->ixblk = pmblk->ixblk;
    pmblk->ixblk = 0;
  }
  if (ipos < 0) {
    m_cache.orgpos = m_cache.mempos;
  } else { /* read from disk */
    ldisk = (long)(ipos - m_cache.nbsav);
    if (ldisk != 0)
      ldisk *= m_cache.lsizb;
    fseek(m_cache.scratch, ldisk, SEEK_SET);
    pdblk = &blk[ipos];
    len = (long)ndel;
    n = (long)pdblk->nsizblk;
    len *= n;
    pvec = pmblk->egyblk;
    lret = 0;
    if (pvec != NULL) {
      lret = (long)fread(pvec, sizeof(double), (size_t)len, m_cache.scratch);
      if (lret == len && (pvec = pmblk->eigblk) != NULL) {
        len = n * n;
        lret = (long)fread(pvec, sizeof(double), (size_t)len, m_cache.scratch);
      }
    }
    if (lret != len)
      throw NumericError("heap read error", CalErrorCode::ScratchFileIo);
    pmblk->nsizblk = pdblk->nsizblk;
    pmblk->ixblk = pdblk->ixblk;
    pdblk->ixblk = 0;
  }
  return pmblk;
}

SBLK *CalCat::sblk_alloc(const int nstruct, const unsigned mxdem)
{
  SBLK *pret, *pblk;
  int k;
  pret = (SBLK *)calalloc((size_t)nstruct * sizeof(SBLK));
  pblk = pret;
  k = nstruct;
  do {
    pblk->eigblk = NULL;
    pblk->egyblk = NULL;
    pblk->ixblk = 0;
    pblk->nsizblk = mxdem;
    ++pblk;
  } while (--k > 0);
  return pret;
}
