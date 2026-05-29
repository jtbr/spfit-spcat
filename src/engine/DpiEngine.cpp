#include "DpiEngine.hpp"
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include "dpi.h"
#include "api/InputSchema.hpp"
#include "common/CalError.hpp"
#include "common/compat.hpp"

DpiEngine::DpiEngine() {
    // Initialize the context with default values
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

int DpiEngine::apply_options(const EngineOptions &opts, int *nfmt, int *itd, int *ndbcd,
                             std::string &namfil)
{
    // DPI option card: "isdgn nvib" (exactly two integers)
    std::ostringstream oss;
    oss << opts.dpi.isdgn << " " << opts.dpi.nvib << "\n";
    std::string card = oss.str();

    FILE *f = fmemopen((void *)card.data(), card.size(), "r");
    if (!f)
        throw IoError("fmemopen failed in DpiEngine::apply_options");

    char namfil_buf[256] = {};
    int n = ::setopt_dpi(&m_context, f, nfmt, itd, ndbcd, namfil_buf);
    fclose(f);

    namfil = namfil_buf;
    return n;
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