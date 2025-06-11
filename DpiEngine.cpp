#include "DpiEngine.hpp"
#include <cstdlib>
#include "dpi.h"

DpiEngine::DpiEngine() {
    // Initialize the context with default values
    m_context.zero = 0.0;
    m_context.nvib = 1; // TODO These initializations seem to be made up, previously unset.
    m_context.iwhole = 0;
    m_context.isdgn = 1;
    m_context.nqn = 4;
}

DpiEngine::~DpiEngine() {
    // No dynamic memory to free
}

int DpiEngine::hamx(int iblk, int nsize, int npar, bcd_t *idpar, double *par,
                      double *egy, double *t, double *dedp,
                      double *pmix, const BOOL ifdump)
{
  return ::hamx_dpi(&m_context, iblk, nsize, npar, idpar, par, egy, t, dedp, pmix, ifdump);
}

int DpiEngine::setint(FILE *lu, BOOL *ifdiag, int *nsav, const int ndip,
                        bcd_t *idip, int *isimag)
{
  return ::setint_dpi(&m_context, lu, ifdiag, nsav, ndip, idip, isimag);
}

int DpiEngine::intens(const int iblk, const int isiz, const int jblk,
                        const int jsiz, const int ndip, const bcd_t *idip,
                        const double *dip, double *s)
{
  return ::intens_dpi(&m_context, iblk, isiz, jblk, jsiz, ndip, idip, dip, s);
}

int DpiEngine::getqn(const int iblk, const int indx, const int maxqn,
                       short *iqn, int *idgn)
{
  return ::getqn_dpi(&m_context, iblk, indx, maxqn, iqn, idgn);
}

int DpiEngine::setopt(FILE *lu, int *nfmt, int *itd, int *ndbcd, char *namfil)
{
  return ::setopt_dpi(&m_context, lu, nfmt, itd, ndbcd, namfil);
}

int DpiEngine::setfmt(int *iqnfmt, int nfmt)
{
  return ::setfmt_dpi(&m_context, iqnfmt, nfmt);
}

int DpiEngine::setblk(FILE *lu, int npar, bcd_t *idpar, double *par,
                        int *nblkpf, int *negy)
{
  return ::setblk_dpi(&m_context, lu, npar, idpar, par, nblkpf, negy);
}