#ifndef _DPI_H_
#define _DPI_H_

// This #ifdef block is the key
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#ifndef BOOL
#define BOOL bool
#endif

/* File contains declarations of public functions whose interface is ~shared with spinv,
   previously in shared header calpgm.h. To be used with context encapsulation by DpiEngine */

typedef unsigned char bcd_t;

/* Function Declarations */
struct DpiContext;

/* Public interface functions, formerly in calpgm.h, now encapsulated in SpinvEngine */
int hamx_dpi(struct DpiContext *ctx, const int iblk, const int nsiz, const int npar, const bcd_t *idpar,
          const double *par, /*@out@*/ double *egy, /*@out@*/ double *t,
          /*@out@*/ double *dedp, /*@out@*/ double *pmix, const BOOL ifdump);
int setint_dpi(struct DpiContext *ctx, FILE *lu, /*@out@*/ BOOL *ifdiag, /*@out@*/ int *nsav, const int ndip,
               bcd_t *idip, /*@out@*/ int *isimag);
int intens_dpi(struct DpiContext *ctx, const int iblk, const int isiz, const int jblk, const int jsiz,
               const int ndip, const bcd_t *idip, const double *dip,
               /*@out@*/ double *s);
int getqn_dpi(struct DpiContext *ctx, const int iblk, const int indx, const int maxqn, /*@out@*/ short *iqn,
              /*@out@*/ int *idgn);
int setopt_dpi(struct DpiContext *ctx, FILE *lu, /*@out@*/ int *nfmt, /*@out@*/ int *itd,
               /*@out@*/ int *nbcdpar, /*@out@*/ char *namfil);
int setfmt_dpi(struct DpiContext *ctx, /*@out@*/ int *iqnfmt, int nfmt);
int setblk_dpi(struct DpiContext *ctx, FILE *lu, const int npar, bcd_t *idpar, const double *par,
               int *nblkpf, int *negy);

#ifdef __cplusplus
} // end extern "C"
#endif

#endif /* _DPI_H_ */