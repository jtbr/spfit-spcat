#ifndef CALCULATION_ENGINE_HPP
#define CALCULATION_ENGINE_HPP

#include <stdio.h>
#include "calpgm.h"

class CalculationEngine {
public:
    virtual ~CalculationEngine() {}

    virtual int hamx(int iblk, int nsize, int npar, bcd_t *idpar, double *par,
                     double *egy, double *teig, double *egyder, double *pmix,
                     int iflg) = 0;
    virtual int setint(int npar, bcd_t *idpar, double *par) = 0;
    virtual double intens(int *iqn, double *dip, int nstat, int *kat,
                          double *egy, double *eigy, int *jrot, int *jroty,
                          int *inv, int *invy, int nsize, int nsiz) = 0;
    virtual void getqn(int iblk, int iqn, int nsize, short *qnum, int *nqn) = 0;
    virtual int setopt(FILE *lu, int *nfmt, int *itd, int *ndbcd, char *namfil) = 0;
    virtual int setfmt(int *iqnfmt, int iflg) = 0;
    virtual int setblk(FILE *lu, int npar, bcd_t *idpar, double *par,
                       int *nblkpf, int *maxdm) = 0;
};

#endif // CALCULATION_ENGINE_HPP
