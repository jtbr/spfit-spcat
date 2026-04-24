/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <climits>
#include "CalCat.hpp"
#include "calpgm.h"
#include "CalError.hpp"
#include "SigintFlag.hpp"

#define PR_DELAY 6
#define NCARD 130

CalCat::CalCat(std::unique_ptr<CalculationEngine> &calc_engine,
               FILE *luout, FILE *lucat, FILE *luegy, FILE *lustr,
               Logger &logger)
    : calc(std::move(calc_engine)),
      m_logger(logger),
      m_luout(luout), m_lucat(lucat), m_luegy(luegy), m_lustr(lustr),
      m_blk(nullptr), m_s(nullptr), m_pmix(nullptr),
      m_par(nullptr), m_derv(nullptr), m_var(nullptr),
      m_idpar(nullptr), m_dip(nullptr), m_idip(nullptr),
      m_nvdip(nullptr), m_isimag(nullptr), m_iqnfmtv(nullptr)
{
  if (!m_luout)
    throw IoError("CalCat constructor received NULL luout stream.", CalErrorCode::FileOpenFailed);
}

CalCat::~CalCat()
{
  if (m_cache.scratch)
    fclose(m_cache.scratch);
}

bool CalCat::run(const CalCatInput &input, CalCatOutput &output)
{
  if (!initializeParameters(input, output))
    return false;
  if (!setupBlocks(input))
    return false;
  if (!computeCatalog(input, output))
    return false;
  if (!finalizeOutput(input, output))
    return false;
  return true;
}

bool CalCat::initializeParameters(const CalCatInput &input, CalCatOutput &output)
{
  m_zero = 1.5e-38;
  m_bigerr = 999.9999;
  m_cmc = 29979.2458;
  m_tmc = 1.43878;

  m_pregy = FALSE;
  m_prstr = FALSE;

  /* Transfer scalars from input */
  m_iflg = input.iflg;
  m_itag = input.itag;
  m_qrot = input.qrot;
  m_inblk = input.inblk;
  m_lblk = input.lblk;
  m_thrsh = input.thrsh;
  m_thrsh1 = input.thrsh1;
  m_fqmax = input.fqmax;
  m_maxv = input.maxv;
  m_ndip = input.ndip;
  m_npdip = input.npdip;
  m_npar = input.npar;
  m_nfit = input.nfit;
  m_nfmt = input.nfmt;
  m_itd = input.itd;
  m_nqn = input.nqn;
  m_catqn = input.catqn;

  /* Build temperature array */
  int ntemp = NTEMP;
  for (int i = 0; i < NTEMP; i++) {
    output.temp[i] = (double)(TMAX - i);
  }
  double tmq = input.tmq;
  for (int k = NTEMP - 2; k >= 0; --k) {
    if (fabs(tmq - output.temp[k]) < 0.01) {
      --ntemp;
      break;
    }
  }
  if (ntemp == NTEMP) {
    output.temp[NTEMP - 1] = tmq;
    for (int k = NTEMP - 2; k >= 0; --k) {
      if (output.temp[k] > tmq)
        break;
      output.temp[k + 1] = output.temp[k];
      output.temp[k] = tmq;
    }
  }
  output.ntemp = ntemp;
  memset(output.qsum, 0, sizeof(output.qsum));
  output.nline = 0;
  output.egymin = 0.;
  output.qrot = m_qrot;

  m_telim = log(1.e-18);

  /* qrot clamping already done in CalCatIO::readInput */

  m_fac = (4.16231e-5) / m_qrot;
  m_prir = (m_iflg >= 1000);
  if (m_prir) {
    m_fac *= m_cmc;
    m_cmc = 1.;
    m_iflg %= 1000;
    m_fqmin = 5e-7;
    m_bigerr *= 0.01;
  } else {
    m_fqmax *= 1000;
    m_fqmin = 5e-5;
  }
  if (m_fqmax < m_fqmin)
    m_fqmax = m_fqmin;

  m_prfrq = (m_iflg >= 100);
  m_iflg %= 100;

  m_prstr = (m_iflg >= 10);
  m_stcomp = 1. / m_zero;
  if (m_prstr) {
    m_stcomp = 1e-7;
    if (m_iflg >= 20)
      m_npdip = -1; /* will be resolved below */
    m_iflg %= 10;
  }
  m_ifdump = (m_iflg >= 5);
  m_prder = (m_iflg == 2 || m_iflg == 4);
  m_preig = (m_iflg >= 3);
  m_pregy = (m_iflg != 0);

  output.ifdump = m_ifdump;

  return true;
}

bool CalCat::setupBlocks(const CalCatInput &input)
{
  size_t nl, nlsq;

  /* Make working copies of input arrays */
  m_par = (double *)calalloc((size_t)m_npar * sizeof(double));
  memcpy(m_par, input.par.data(), (size_t)m_npar * sizeof(double));

  m_derv = (double *)calalloc((size_t)m_npar * sizeof(double));
  memcpy(m_derv, input.derv.data(), (size_t)m_npar * sizeof(double));

  nl = (size_t)(m_npar * input.ndbcd);
  m_idpar = (bcd_t *)calalloc(nl);
  memcpy(m_idpar, input.idpar.data(), nl);

  m_var = (double *)calalloc(input.var.size() * sizeof(double));
  memcpy(m_var, input.var.data(), input.var.size() * sizeof(double));

  m_dip = (double *)calalloc((size_t)m_ndip * sizeof(double));
  memcpy(m_dip, input.dip.data(), (size_t)m_ndip * sizeof(double));

  nl = (size_t)(m_ndip * NDECDIP);
  m_idip = (bcd_t *)calalloc(nl);
  memcpy(m_idip, input.idip.data(), nl);

  m_nvdip = (int *)calalloc((size_t)m_ndip * sizeof(int));
  memcpy(m_nvdip, input.nvdip.data(), (size_t)m_ndip * sizeof(int));

  m_isimag = (int *)calalloc((size_t)m_ndip * sizeof(int));
  memcpy(m_isimag, input.isimag.data(), (size_t)m_ndip * sizeof(int));

  m_iqnfmtv = (int *)calalloc((size_t)(m_nfmt << 1) * sizeof(int));
  memcpy(m_iqnfmtv, input.iqnfmtv.data(), (size_t)(m_nfmt << 1) * sizeof(int));

  /* Resolve npdip from input — re-derive from npdip and prstr flags */
  /* The npdip logic was partially handled in CalCatIO, but the prstr-based
     override happens here since prstr depends on iflg which we already parsed */
  if (m_prstr) {
    if (m_npdip < 0) {
      m_npdip = -m_npdip; /* was set negative by CalCatIO to signal multi-component */
    }
  } else {
    m_nvdip[0] = m_ndip;
    m_npdip = 1;
  }

  m_ndel = (unsigned)(m_nfit + 1);

  /* QN format processing */
  m_iqnfmt = m_iqnfmtv[0];
  m_iposv = 0;
  if (m_nfmt > 1 || ((m_iqnfmt / 1000) & 1) != 0)
    m_iposv = ((m_iqnfmt / 100) % 5) - 1;
  if (m_iposv == 0 || m_maxv < 0)
    m_maxv = 1000;
  m_globfmt = 0;
  for (int i = 0; i < m_nfmt; ++i) {
    int j = m_iqnfmtv[i] / 1000;
    if (j != 4 && j != 5)
      j = 0;
    m_iqnfmtv[i + m_nfmt] = j;
    m_globfmt += j;
  }
  m_newfmt = m_iqnfmtv[m_nfmt];

  int nsize_p = INT_MAX;
  if (m_nfit > nsize_p || m_nfit <= 0) {
    m_logger.error("var matrix size invalid (nfit = %d).", m_nfit);
    return false;
  }

  m_nbkpj = m_lblk;
  int i_maxdm = nsize_p;
  fprintf(m_luout, "IQNFMT = %d\n", m_iqnfmt);
  calc->setblk(m_luout, m_npar, m_idpar, m_par, &m_nbkpj, &i_maxdm);
  m_maxdm = (unsigned)i_maxdm;

  m_inblk = m_nbkpj * m_inblk + 1;
  m_lblk = m_nbkpj * (m_lblk + 1);

  nl = (size_t)m_maxdm * sizeof(double);
  m_pmix = (double *)calalloc(nl);
  m_pmix[0] = 1.;

  nlsq = nl * (size_t)m_maxdm;
  m_s = (double **)calalloc((size_t)m_npdip * sizeof(double *));
  m_s[0] = (double *)calalloc(nlsq);
  for (int k = 1; k < m_npdip; ++k) {
    m_s[k] = (double *)calalloc(nlsq);
  }

  /* set up DIAG, NSAV */
  calc->setint(m_luout, &m_diag, &m_nsav, m_ndip, m_idip, m_isimag);

  if (m_prstr) {
    int k = 0, ij = 0;
    int i_end = m_nvdip[0];
    for (int j = 1; j < m_ndip; ++j) {
      if (j == i_end) {
        m_isimag[++k] = m_isimag[j];
        i_end += m_nvdip[k];
      } else if (m_isimag[j] != m_isimag[k]) {
        if (m_isimag[k] < 0) {
          m_isimag[k] = m_isimag[j];
        } else if (m_isimag[j] >= 0) {
          ++ij;
        }
      }
    }
    if (ij > 0)
      fputs(" mixed magnetic and electric dipoles in str file\n", m_luout);
  }

  m_nsav = (m_nsav + 1) * m_nbkpj;
  if (m_diag && (m_ifdump || m_maxdm == 1))
    m_diag = 0;
  if (!m_ifdump && m_preig)
    m_preig = m_diag;

  /* allocate block structure */
  m_blk = sblk_alloc(m_nsav + 1, m_maxdm);
  nl = (size_t)m_maxdm;
  nl = ((size_t)m_ndel + nl) * nl;
  {
    /* Assume plenty of address space; cap only against NDHEAPC if defined. */
    m_nbsav = m_nsav;
  }
  nlsq = nl * sizeof(double);
#ifdef NDHEAPC
#if    NDHEAPC
  {
    size_t ndh = (size_t)NDHEAPC;
    ndh = ndh / nlsq;
    if (ndh < (size_t)m_nbsav)
      m_nbsav = (int)ndh;
  }
#endif
#endif
  if (m_nbsav < 2)
    m_nbsav = 2;

  nl = (size_t)m_maxdm * sizeof(double);
  nlsq = nl * (size_t)m_maxdm;
  nl *= (size_t)m_ndel;

  SBLK *pblk = m_blk;
  for (int i = 0; i <= m_nbsav; ++i) {
    pblk->egyblk = (double *)calalloc(nl);
    if (m_diag) {
      pblk->eigblk = (double *)calalloc(nlsq);
    }
    ++pblk;
  }
  pblk = NULL;
  /* set up intensity constants */
  m_tmc = -m_tmc;
  m_tmq = m_tmc / (input.tmq * m_cmc);
  m_tmql = m_tmq * 0.43429448; /* tmql = log10(exp(1))*tmq */
  m_thrsh1 -= 10.9542425;
  m_starg = m_thrsh - log10(m_fqmax * m_fac);
  double scomp = -38;
  m_strmn = (m_starg > scomp) ? m_starg : scomp;
  m_strmn = pow(10., m_strmn);

  short iqni[MAXQN + MAXQNX];
  int isiz;
  m_maxqn = calc->getqn(m_inblk, 0, 0, iqni, &isiz);
  if (m_maxqn > MAXQNX)
    m_maxqn = MAXQNX;

  return true;
}

bool CalCat::computeCatalog(const CalCatInput &/*input*/, CalCatOutput &output)
{
  char sqn[4 * MAXQN + 2];
  memset(sqn, 0, sizeof(sqn));
  char sgup[4];
  short iqni[MAXQN + MAXQNX];
  int isiz, jsiz, idgn, jdgn, kbgn;
  double *teig, *teigp, *egy, *egyp, *dedp, *sij;
  double dgn, dgnf, egx, err, frq, str, strr, strq, strlg, elow, te, tmp;
  double diff, thrshf;
  long nline = 0;
  double egymin = 0.;
  int ntemp = output.ntemp;
  int ktsp = -1;
  BOOL first;
  int igup, isneg;

  int catqn;
  if (m_nqn > MAXCAT) {
    catqn = m_nqn;
  } else {
    sqn[24] = ' ';
    catqn = m_catqn;
  }
  catqn = catqn << 1;
  memset(sqn, (int)' ', (size_t)(catqn << 1));

  /* START MAJOR LOOP OVER BLOCKS */
  for (int iblk = m_inblk; iblk <= m_lblk; ++iblk) {
    first = m_prfrq;
    int j = (iblk - 1) / m_nbkpj;
    if ((iblk - j * m_nbkpj) == 1) {
      if (SigintFlag::isTriggered())
        break;
      SBLK *pblk = m_blk;
      int jblk = iblk + m_nbkpj - m_nsav;
      for (int i = 0; i <= m_nsav; ++i) {
        if (pblk->ixblk < jblk)
          pblk->ixblk = 0;
        ++pblk;
      }
      if (caldelay(PR_DELAY)) {
        m_logger.info(" STARTING QUANTUM %3d", j);
        fflush(stdout);
      }
    }
    /* move first block to unused position */
    SBLK *pblk = ibufof(-1, m_ndel, m_blk);
    /* get energy and derivatives for new block */
    calc->getqn(iblk, 0, 0, iqni, &isiz);
    if (isiz <= 0)
      continue;
    if ((unsigned)isiz > m_maxdm)
      break;
    pblk->ixblk = iblk;
    pblk->nsizblk = (unsigned)isiz;
    teig = NULL;
    if (m_diag)
      teig = pblk->eigblk;
    if (teig == NULL)
      teig = m_s[0];
    egy = pblk->egyblk;
    if (egy == NULL)
      break;
    dedp = egy + isiz;
    calc->hamx(iblk, isiz, m_npar, m_idpar, m_par, egy, teig, dedp, m_pmix, m_ifdump);

    /* print out energies and compute partition function */
    if (m_prfrq) {
      fprintf(m_luout, " ENERGIES FOR BLOCK NUMBER %3d, INDEX-DEGEN-ENERGY-",
              iblk);
      fputs("-EST.ERROR-MIXING-   QUANTUM NUMBERS\n", m_luout);
      first = TRUE;
    }
    kbgn = 0;
    teigp = teig;
    for (int i = 0; i < isiz; ++i, teigp += isiz) {
      if (m_ifdump)
        kbgn = i;
      calc->getqn(iblk, i + 1, m_nqn, iqni, &idgn);
      if (idgn <= 0)
        continue;
      int iv = iqni[m_iposv];
      if (iv > m_maxv)
        continue;
      if (i == 0 && m_nfmt > 1 && iv < m_nfmt) {
        m_iqnfmt = m_iqnfmtv[iv];
        m_newfmt = m_iqnfmtv[iv + m_nfmt];
      }
      if (m_globfmt != 0) {
        ktsp = -1;
        if (m_newfmt != 0)
          ktsp = iqni[m_nqn - 2];
        calc->getqn(iblk, i + 1, m_maxqn, iqni, &idgn);
      }
      if (m_diag) {
        if (fabs(teigp[i]) < 0.01) {
          int k = (int)idamax(isiz, teigp, 1);
          if (fabs(teigp[k]) < 0.01)
            idgn = 0;
        }
      }
      /* calculate uncertainty */
      dcopy(m_nfit, &dedp[i], isiz, m_derv, 1);
      err = calerr(m_nfit, m_var, m_derv) / m_cmc;
      if (err > m_bigerr)
        err = m_bigerr;
      egx = egy[i] / m_cmc;
      if (egymin > egx) { /* adjust energy zero */
        te = m_tmc * (egymin - egx);
        egymin = egx;
        for (int k = 0; k < ntemp; ++k) {
          tmp = te / output.temp[k];
          if (tmp < m_telim)
            tmp = m_telim;
          output.qsum[k] *= exp(tmp);
        }
      }
      /* accumulate partition function values */
      te = m_tmc * (egx - egymin);
      dgn = (double)idgn;
      for (int k = 0; k < ntemp; ++k) {
        tmp = te / output.temp[k];
        if (tmp < m_telim)
          break;
        output.qsum[k] += dgn * exp(tmp);
      }
      if (m_prfrq) {
        fprintf(m_luout, "%6d %4d %17.6f %17.6f %10.6f",
                iblk, i + 1, egx, err, m_pmix[i]);
        if (m_globfmt != 0)
          fprintf(m_luout, " (%3d)", ktsp);
        for (int k = 0; k < m_maxqn; ++k) {
          fprintf(m_luout, "%3d", (int)iqni[k]);
        }
        fputc('\n', m_luout);
      }
      if (m_pregy) {
        fprintf(m_luegy, "%6d %4d %17.6f %17.6f %10.6f %4d:",
                iblk, i + 1, egx, err, m_pmix[i], idgn);
        if (m_globfmt != 0)
          fprintf(m_luegy, "(%3d)", ktsp);
        for (int k = 0; k < m_maxqn; ++k) {
          fprintf(m_luegy, "%3d", (int)iqni[k]);
        }
        fputc('\n', m_luegy);
      }
      if (m_prder) {
        j = 0;
        for (int k = 0; k < m_nfit; ++k) {
          if (j == 0)
            fputs("        ", m_luegy);
          fprintf(m_luegy, " %5d %13.5E", k + 1, m_derv[k]);
          if ((++j) == 6) {
            fputc('\n', m_luegy);
            j = 0;
          }
        }
        if (j > 0)
          fputc('\n', m_luegy);
      }
      if (m_preig) {
        j = 0;
        for (int k = kbgn; k < isiz; ++k) {
          if (j == 0)
            fputs("        ", m_luegy);
          fprintf(m_luegy, " %5d %13.5E", k + 1, teigp[k]);
          if ((++j) == 6) {
            fputc('\n', m_luegy);
            j = 0;
          }
        }
        if (j > 0)
          fputc('\n', m_luegy);
      }
    }

    /* loop over previous blocks for intensity */
    int jblk = 0;
    do {
      int jnxt = iblk + 1;
      pblk = m_blk;
      int jj = 0;
      for (j = 0; j <= m_nsav; ++j) {
        int jx = pblk->ixblk;
        if (jx > jblk && jx < jnxt) {
          jnxt = jx;
          jj = j;
        }
        ++pblk;
      }
      jblk = jnxt;
      pblk = ibufof(jj, m_ndel, m_blk);
      jsiz = (int)pblk->nsizblk;
      teigp = pblk->eigblk;
      if (teigp == NULL)
        teigp = m_s[0];
      idgn = 0;
      int ij = 0;
      for (int i = 0; i < m_npdip; ++i) {
        sij = m_s[i];
        int idf = calc->intens(iblk, isiz, jblk, jsiz, m_nvdip[i],
                               &m_idip[ij * NDECDIP], &m_dip[ij], sij);
        if (idf != 0) {
          if (m_diag)
            simt(isiz, jsiz, sij, teig, teigp, m_pmix);
          idgn = idf;
        }
        sij = NULL;
        ij += m_nvdip[i];
      }
      if (idgn <= 0)
        continue;
      dgnf = (double)idgn;
      egyp = pblk->egyblk;
      if (egyp == NULL)
        continue;
      int idf = (int)m_maxdm;
      int jmax = jsiz;
      /* Need a dvec for str output */
      double dvec_local[10]; /* enough for npdip <= 10 */
      double *dvec_str = dvec_local;
      double *dvec_alloc = NULL;
      if (m_npdip > 10) {
        dvec_alloc = (double *)calalloc((size_t)m_npdip * sizeof(double));
        dvec_str = dvec_alloc;
      }
      for (int i = 0; i < isiz; ++i) {
        calc->getqn(iblk, i + 1, m_nqn, iqni, &idgn);
        if (idgn <= 0 || iqni[m_iposv] > (short)m_maxv)
          continue;
        dgn = idgn / dgnf;
        if (iblk == jblk)
          jmax = i + 1;
        for (j = 0; j < jmax; ++j) {
          if (i - j + idf < 0)
            break;
          calc->getqn(jblk, j + 1, m_nqn, &iqni[MAXQN], &jdgn);
          if (jdgn <= 0 || iqni[m_iposv + MAXQN] > (short)m_maxv)
            continue;
          frq = egy[i] - egyp[j];
          int ixdp = i;
          int jxdp = j;
          for (int k = 0; k < m_nfit; ++k) {
            ixdp += isiz;
            jxdp += jsiz;
            m_derv[k] = egy[ixdp] - egyp[jxdp];
          }
          isneg = 0;
          if (frq < 0.) {
            isneg = 1;
            frq = -frq;
            qnfmt(&iqni[MAXQN], m_nqn, sqn);
            qnfmt(iqni, m_nqn, sqn + catqn);
            igup = jdgn;
            elow = egy[i];
          } else {
            qnfmt(iqni, m_nqn, sqn);
            qnfmt(&iqni[MAXQN], m_nqn, sqn + catqn);
            igup = idgn;
            elow = egy[i] - frq;
          }
          long ioff = isiz;
          ioff = i + ioff * j;
          strr = 0.;
          strq = 0.;
          for (ij = 0; ij < m_npdip; ++ij) {
            sij = m_s[ij];
            str = sij[ioff];
            dvec_str[ij] = str;
            strr += str;
            strq += str * str;
          }
          if (strq > m_stcomp) {
            for (ij = 0; ij < m_npdip; ++ij) {
              str = dvec_str[ij];
              if (isneg != 0 && m_isimag[ij] > 0)
                str = -str;
              fprintf(m_lustr, "%15.4f%15.6E %4d %s %4d\n", frq, str,
                      m_iqnfmt, sqn, ij + 1);
            }
          }
          strr *= strr;
          if (frq > m_fqmax)
            continue;
          if (frq < m_fqmin)
            continue;
          if (strr < m_strmn)
            continue;
          str = dgn * strr * m_fac * frq * (1. - exp(m_tmq * frq));
          strlg = log10(str + m_zero) + m_tmql * elow;
          if (strlg < m_thrsh)
            continue;
          thrshf = m_thrsh1 + log10(frq + m_zero) * 2.;
          diff = m_thrsh - thrshf;
          if (fabs(diff) < 4.) {
            diff = pow(10., diff) + 1.;
            thrshf = log10(diff) + thrshf;
          }
          if (strlg < thrshf)
            continue;
          /* calculate errors */
          err = calerr(m_nfit, m_var, m_derv);
          if (err > m_bigerr)
            err = m_bigerr;
          elow /= m_cmc;
          if (first) {
            fputs(" FREQUENCY-EST.ERROR.-LINE.STR. DIP**2-LGSTR.-ITD,",
                  m_luout);
            fputs("-GUP-I.D.-QNFORM-QUANTUM NUMBERS\n", m_luout);
            first = FALSE;
          }
          if (m_prfrq) {
            fprintf(m_luout, "%13.4f %8.4f %12.5E %8.4f %2d",
                    frq, err, strr, strlg, m_itd);
            fprintf(m_luout, "%10.4f %3d %7ld %4d %s\n",
                    elow, igup, m_itag, m_iqnfmt, sqn);
          }
          ++nline;
          gupfmt(igup, sgup);
          sgup[3] = '\0';
          if (m_prir) {
            fprintf(m_lucat, "%13.6f%8.6f%8.4f%2d%10.4f%s%7ld%4d",
                    frq, err, strlg, m_itd, elow, sgup, m_itag, m_iqnfmt);
          } else if (frq < 99999999.) {
            fprintf(m_lucat, "%13.4f%8.4f%8.4f%2d%10.4f%s%7ld%4d",
                    frq, err, strlg, m_itd, elow, sgup, m_itag, m_iqnfmt);
          } else {
            fprintf(m_lucat, "%13.3f%8.3f%8.3f%2d%10.4f%s%7ld%4d",
                    frq, err, strlg, m_itd, elow, sgup, m_itag, m_iqnfmt);
          }
          fputs(sqn, m_lucat);
          fputc('\n', m_lucat);
        }
      }
      if (dvec_alloc)
        free(dvec_alloc);
    } while (jblk != iblk);
  }
  teig = teigp = NULL;

  output.nline = nline;
  output.egymin = egymin;

  return true;
}

bool CalCat::finalizeOutput(const CalCatInput &/*input*/, CalCatOutput &output)
{
  static char warn[] = " WARNING: THERE WAS NO DIAGONALIZATION\n";
  static char headq[] = "TEMPERATURE - Q(SPIN-ROT.) - log Q(SPIN-ROT.)\n";

  if (m_ifdump) {
    m_logger.warn(" WARNING: THERE WAS NO DIAGONALIZATION");
    fputs(warn, m_luout);
  }
  m_logger.info("INITIAL Q = %14.4f, NEW Q IS RELATIVE TO MIN.EGY.= %14.4f",
                m_qrot, output.egymin);
  fprintf(m_luout, "INITIAL Q = %14.4f, NEW Q IS RELATIVE TO MIN.EGY.= %14.4f\n",
          m_qrot, output.egymin);
  m_logger.info(" NUMBER OF LINES = %6ld", output.nline);
  fprintf(m_luout, " NUMBER OF LINES = %6ld\n", output.nline);
  m_logger.info("TEMPERATURE - Q(SPIN-ROT.) - log Q(SPIN-ROT.)");
  fputs(headq, m_luout);
  for (int i = 0; i < output.ntemp; ++i) {
    double qlog = -100;
    if (output.qsum[i] > m_zero)
      qlog = log10(output.qsum[i]);
    m_logger.info(" %10.3f %14.4f %9.4f", output.temp[i], output.qsum[i], qlog);
    fprintf(m_luout, " %10.3f %14.4f %9.4f\n", output.temp[i], output.qsum[i], qlog);
  }

  /* release engine storage */
  int nbkpj_tmp = m_nbkpj;
  int i_tmp = 0;
  calc->setblk(m_luout, 0, m_idpar, m_par, &nbkpj_tmp, &i_tmp);

  /* free block structure memory */
  SBLK *pblk = m_blk;
  for (int i = 0; i <= m_nsav; ++i) {
    if (pblk->eigblk != NULL)
      free(pblk->eigblk);
    if (pblk->egyblk != NULL)
      free(pblk->egyblk);
    ++pblk;
  }
  free(m_blk);
  m_blk = NULL;

  /* free intensity matrices */
  for (int i = 0; i < m_npdip; ++i) {
    free(m_s[i]);
    m_s[i] = NULL;
  }
  free(m_s);
  m_s = NULL;
  free(m_pmix);
  m_pmix = NULL;

  /* free working copies */
  if (m_var != NULL)
    free(m_var);
  free(m_idpar);
  free(m_derv);
  free(m_par);
  free(m_iqnfmtv);
  free(m_idip);
  free(m_isimag);
  free(m_nvdip);
  free(m_dip);

  m_var = m_par = m_derv = m_dip = m_pmix = NULL;
  m_idpar = m_idip = NULL;
  m_nvdip = m_isimag = m_iqnfmtv = NULL;

  return true;
}
