#include <cstring>
#include <cstdlib>
#include "spinv_internal.h"
#include "SpinvEngine.hpp"

SpinvEngine::SpinvEngine()
{
  // Zero out the entire context first
  memset(&m_context, 0, sizeof(SpinvContext));

  // Apply specific initializations
  m_context.sptzero[0] = 1;
  m_context.sptzero[1] = 0;
  m_context.ssp_head.next = NULL; // ssp_head itself needs init
  m_context.ssp_head.sspt = m_context.sptzero;
  m_context.ssp_head.ssize = 1;
  m_context.ssp_head.nitot = 0;
  m_context.revsym[0] = 0;
  m_context.revsym[1] = 3;
  m_context.revsym[2] = 2;
  m_context.revsym[3] = 1;
  m_context.isoddk[0] = 0;
  m_context.isoddk[1] = 1;
  m_context.isoddk[2] = 1;
  m_context.isoddk[3] = 0;
  m_context.ixphase[0] = 0;
  m_context.ixphase[1] = 1;
  m_context.ixphase[2] = 2;
  m_context.ixphase[3] = 3;
  m_context.ipwr2[0] = 0;
  m_context.ipwr2[1] = 1;
  m_context.ipwr2[2] = 2;
  m_context.ipwr2[3] = 4;

  // Initialize pointers to static fallback members
  m_context.vinfo = &m_context.vinfo1;        // vinfo1 itself needs full init
  memset(&m_context.vinfo1, 0, sizeof(SVIB)); // Zero vinfo1
  m_context.vinfo1.spt = m_context.sptzero;   // Then set its specific pointers

  m_context.dipinfo = &m_context.dipinfo0; // dipinfo0 needs full init
  memset(&m_context.dipinfo0, 0, sizeof(SDIP));

  // Pointing arrays to their single static "zero" element is okay for initial state
  // as setopt/setblk will realloc them.
  m_context.moldv = &m_context.zmoldv;
  m_context.blkptr = &m_context.zblkptr;
  m_context.ipder = &m_context.zipder;
  m_context.ivs = &m_context.zivs;
  m_context.wk = &m_context.zwk;
  m_context.idx = &m_context.zidx;
  m_context.jdx = &m_context.zjdx;
  m_context.iqnsep = &m_context.ziqnsep;
  m_context.ibkptr = &m_context.zibkptr;
  m_context.ikmin = &m_context.zikmin;

  // The single static elements themselves
  m_context.zwk = 0.0;
  m_context.zmoldv = 0;
  m_context.zblkptr = 0;
  m_context.zivs = 0;
  m_context.zipder = 0; // ipder values are significant; -1 is often "not fitted"
  m_context.zidx = 0;
  m_context.zjdx = 0;
  m_context.ziqnsep = 0;
  m_context.zibkptr = 0;
  m_context.zikmin = 0;

  // Initialize glob explicitly too (setopt overwrites many but not all)
  memset(&m_context.glob, 0, sizeof(GLOB));
  // setopt then sets: glob.lsym=TRUE, glob.esym=TRUE, glob.esymdec=100 etc.
}

SpinvEngine::~SpinvEngine()
{
  // Free any dynamically allocated memory in the context
  if (m_context.vinfo != &m_context.vinfo1) {
      free(m_context.vinfo);
  }
  if (m_context.moldv != &m_context.zmoldv) {
      free(m_context.moldv);
  }
  if (m_context.blkptr != &m_context.zblkptr) {
      free(m_context.blkptr);
  }
  if (m_context.ipder != &m_context.zipder) {
      free(m_context.ipder);
  }
  if (m_context.ivs != &m_context.zivs) {
      free(m_context.ivs);
  }
  if (m_context.wk != &m_context.zwk) {
      free(m_context.wk);
  }
  if (m_context.idx != &m_context.zidx) {
      free(m_context.idx);
  }
  if (m_context.jdx != &m_context.zjdx) {
      free(m_context.jdx);
  }
  if (m_context.iqnsep != &m_context.ziqnsep) {
      free(m_context.iqnsep);
  }
  if (m_context.ibkptr != &m_context.zibkptr) {
      free(m_context.ibkptr);
  }
  if (m_context.ikmin != &m_context.zikmin) {
      free(m_context.ikmin);
  }
}

int SpinvEngine::hamx(int iblk, int nsize, int npar, bcd_t *idpar, double *par,
                      double *egy, double *t, double *dedp,
                      double *pmix, const BOOL ifdump)
{
    return ::hamx(&m_context, iblk, nsize, npar, idpar, par, egy, t, dedp, pmix, ifdump);
}

int SpinvEngine::setint(FILE *lu, BOOL *ifdiag, int *nsav, const int ndip,
                        bcd_t *idip, int *isimag)
{
    return ::setint(&m_context, lu, ifdiag, nsav, ndip, idip, isimag);
}

int SpinvEngine::intens(const int iblk, const int isiz, const int jblk,
                           const int jsiz, const int ndip, const bcd_t *idip,
                           const double *dip, double *s)
{
    return ::intens(&m_context, iblk, isiz, jblk, jsiz, ndip, idip, dip, s);
}

int SpinvEngine::getqn(const int iblk, const int indx, const int maxqn,
                        short *iqn, int *idgn)
{
    return ::getqn(&m_context, iblk, indx, maxqn, iqn, idgn);
}

int SpinvEngine::setopt(FILE *lu, int *nfmt, int *itd, int *ndbcd, char *namfil)
{
    printf("DEBUG setopt at entry: m_context.glob.nqn0 = %d\n", m_context.glob.nqn0);
    int retval = ::setopt(&m_context, lu, nfmt, itd, ndbcd, namfil);
    printf("DEBUG setopt at exit: m_context.glob.nqn0 = %d\n", m_context.glob.nqn0);
    return retval;
}

int SpinvEngine::setfmt(int *iqnfmt, int nfmt)
{
    return ::setfmt(&m_context, iqnfmt, nfmt);
}

int SpinvEngine::setblk(FILE *lu, int npar, bcd_t *idpar, double *par,
                        int *nblkpf, int *negy)
{
    return ::setblk(&m_context, lu, npar, idpar, par, nblkpf, negy);
}