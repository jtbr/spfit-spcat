#include "SpinvEngine.hpp"
#include "spinv_internal.h"
#include <cstdlib>

SpinvEngine::SpinvEngine() :
  m_context {
    .sptzero = {1, 0},
    .ssp_head = {NULL, NULL, 1, 0},
    .revsym = {0, 3, 2, 1},
    .isoddk = {0, 1, 1, 0},
    .ixphase = {0, 1, 2, 3},
    .ipwr2 = {0, 1, 2, 4}
  }
{
    // Initialize the context with default values
    // Pointers
    m_context.vinfo = &m_context.vinfo1;
    m_context.ssp_head.sspt = m_context.sptzero;
    m_context.dipinfo = &m_context.dipinfo0;
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

    // initial values
    m_context.zero = 0.;
    m_context.zwk = 0.;
    m_context.zmoldv = 0;
    m_context.zblkptr = 0;
    m_context.zivs = 0;
    m_context.zipder = 0;
    m_context.szero = 0;
    m_context.zidx = 0;
    m_context.zjdx = 0;
    m_context.ziqnsep = 0;
    m_context.zibkptr = 0;
    m_context.zikmin = 0;
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
                      double *egy, double *teig, double *egyder, double *pmix,
                      int iflg) {
    // This function will be implemented later, once the global variables
    // in the spinv_*.c files have been replaced with the context struct.
    return 0;
}

int SpinvEngine::setint(int npar, bcd_t *idpar, double *par) {
    // This function will be implemented later.
    return 0;
}

double SpinvEngine::intens(int *iqn, double *dip, int nstat, int *kat,
                           double *egy, double *eigy, int *jrot, int *jroty,
                           int *inv, int *invy, int nsize, int nsiz) {
    // This function will be implemented later.
    return 0.0;
}

void SpinvEngine::getqn(int iblk, int iqn, int nsize, short *qnum, int *nqn) {
    // This function will be implemented later.
}

int SpinvEngine::setopt(FILE *lu, int *nfmt, int *itd, int *ndbcd, char *namfil) {
    // This function will be implemented later.
    return 0;
}

int SpinvEngine::setfmt(int *iqnfmt, int iflg) {
    // This function will be implemented later.
    return 0;
}

int SpinvEngine::setblk(FILE *lu, int npar, bcd_t *idpar, double *par,
                        int *nblkpf, int *maxdm) {
    // This function will be implemented later.
    return 0;
}
