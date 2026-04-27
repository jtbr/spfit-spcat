#ifndef _ULIB_H_
#define _ULIB_H_

#include "calpgm_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/************** ULIB interface ******************************************/
int ordblk(const int ndm, const int nn, short *iqnsep, double *t, double *e,
            short *isblk, double *p, short *ip);
void etswap(const int ndm,const int nsize,const int ix1,const int ix2,
            double *t,/*@null@*/double *e, /*@null@*/double *q);
int hdiag(const int nm, const int nx, double *z, /*@out@*/ double *d,
          /*@out@*/ double *e, /*@out@*/ short *iqnsep);
int triag(const int nm, const int nx, const int nz,
          double *z, double *d, double *e);
int getpar(FILE *lu, FILE *luout, /*@out@*/int *npar, /*@out@*/int* npartot,
           /*@out@*/ bcd_t *idpar,/*@out@*/double *par,
           /*@out@*/ double *erpar, /*@out@*/char *plbl, int plblen);
int getvar(FILE *lu, const int npar, /*@out@*/ double *var, bcd_t *idpar,
           double *erpar, const int flg);
int putvar(FILE *lu, const int npar, double *var, double * erpar);
double calerr(const int npar,const double *var,const double *derv);
int deflin(int iqnfmt, /*@out@*/ short *idqn);
int getlin(FILE *lu, const int nqn, short *idqn, /*@out@*/ short *iqn,
           /*@out@*/ double *xfreq, /*@out@*/ double *xerr,
           /*@out@*/ double *xwt, /*@out@*/ char *card, const int ncard);
int getbcd(const char *line, bcd_t *ivbcd, int nbcd);
int putbcd(/*@out@*/ char *line, int nlen, const bcd_t *ivbcd);
int bcd2i(bcd_t btmp);
bcd_t i2bcd(int i);

#ifdef __cplusplus
}
#endif

#endif /* _ULIB_H_ */
