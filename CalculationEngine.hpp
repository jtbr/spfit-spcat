#ifndef CALCULATION_ENGINE_HPP
#define CALCULATION_ENGINE_HPP

#include <stdio.h>
#include "calpgm.h"

class CalculationEngine {
public:
    virtual ~CalculationEngine() {}

    virtual int hamx(int iblk, int nsize, int npar, bcd_t *idpar, double *par,
                     double *egy, double *t,
                     double *dedp, double *pmix, const BOOL ifdump) = 0;
    virtual int setint(FILE *lu, BOOL *ifdiag, int *nsav, const int ndip,
                       bcd_t *idip, int *isimag) = 0;
    virtual int intens(const int iblk, const int isiz, const int jblk,
                       const int jsiz, const int ndip, const bcd_t *idip,
                       const double *dip, double *s) = 0;
    virtual int getqn(const int iblk, const int indx, const int maxqn,
                      short *iqn, int *idgn) = 0;
    virtual int setopt(FILE *lu, int *nfmt, int *itd, int *ndbcd, char *namfil) = 0;
    virtual int setfmt(int *iqnfmt, int nfmt) = 0;
    virtual int setblk(FILE *lu, int npar, bcd_t *idpar, double *par,
                       int *nblkpf, int *negy) = 0;
};

#endif // CALCULATION_ENGINE_HPP
