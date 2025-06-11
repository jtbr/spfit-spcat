#ifndef DPI_ENGINE_HPP
#define DPI_ENGINE_HPP

#include "CalculationEngine.hpp"
#include "DpiContext.hpp"

class DpiEngine : public CalculationEngine {
  public:
    DpiEngine();
    ~DpiEngine();

    int hamx(int iblk, int nsize, int npar, bcd_t *idpar, double *par,
             double *egy, double *t, double *dedp,
             double *pmix, const BOOL ifdump) override;
    int setint(FILE *lu, BOOL *ifdiag, int *nsav, const int ndip,
               bcd_t *idip, int *isimag) override;
    int intens(const int iblk, const int isiz, const int jblk,
                  const int jsiz, const int ndip, const bcd_t *idip,
                  const double *dip, double *s) override;
    int getqn(const int iblk, const int indx, const int maxqn,
               short *iqn, int *idgn) override;
    int setopt(FILE *lu, int *nfmt, int *itd, int *ndbcd, char *namfil) override;
    int setfmt(int *iqnfmt, int nfmt) override;
    int setblk(FILE *lu, int npar, bcd_t *idpar, double *par,
               int *nblkpf, int *negy) override;

  private:
    DpiContext m_context;
};

#endif // DPI_ENGINE_HPP
