#ifndef DPI_ENGINE_HPP
#define DPI_ENGINE_HPP

#include "CalculationEngine.hpp"
#include "DpiContext.hpp"

class DpiEngine : public CalculationEngine {
public:
    DpiEngine();
    ~DpiEngine();

    int hamx(int iblk, int nsize, int npar, bcd_t *idpar, double *par,
             double *egy, double *teig, double *egyder, double *pmix,
             int iflg) override;
    int setint(int npar, bcd_t *idpar, double *par) override;
    double intens(int *iqn, double *dip, int nstat, int *kat,
                  double *egy, double *eigy, int *jrot, int *jroty,
                  int *inv, int *invy, int nsize, int nsiz) override;
    void getqn(int iblk, int iqn, int nsize, short *qnum, int *nqn) override;
    int setopt(FILE *lu, int *nfmt, int *itd, int *ndbcd, char *namfil) override;
    int setfmt(int *iqnfmt, int iflg) override;
    int setblk(FILE *lu, int npar, bcd_t *idpar, double *par,
               int *nblkpf, int *maxdm) override;

private:
    DpiContext m_context;
};

#endif // DPI_ENGINE_HPP
