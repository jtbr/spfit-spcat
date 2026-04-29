#ifndef _SUBFIT_H_
#define _SUBFIT_H_

#include "splib/calpgm_types.h"

/************** SUBFIT interfaces ***********************************/
/* subfit.cpp is C++; all callers are C++, no extern "C" needed. */
int getlbl(int npar, bcd_t *idpar, /*@out@*/ char *parlbl, const char *fil, int idiv,
           int len);
int filbak(char *flu, char *fbak);
int prcorr(FILE *lufit, int npar, double *cor, int ndcor, double *err,
           int normalize_to_correlation);
SXLINE *line_at(int ipos);
void init_line_buffer(size_t nlines, int nfit);
void release_line_buffer();
void dnuadd(int npar, int nparx, int initl, int indx,
            int ifac, double *egy, double *egyder, int nsize,
            int line, const double *par, const double *fac);
double dnuget(int iflg, int npar, double f, int line, double *dvec);
int getdbk(int *link, /*@out@*/ int *iblk, /*@out@*/ int *indx,
           /*@out@*/ int *initl, /*@out@*/ int *ifac);
int frqdat(int line, /*@out@*/ int *ibln, /*@out@*/ double *xfrq,
        /*@out@*/ double *xwt, /*@out@*/ double *xerr, /*@out@*/ short *iqn);
int lnlink(int *prvblk, int nblk, int iblk, int line);

#endif /* _SUBFIT_H_ */
