/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
#ifndef _CALPGM_TYPES_H_
#define _CALPGM_TYPES_H_

#include <stdio.h>
#include <stdbool.h>
#ifndef BOOL
#define BOOL bool
#endif
#ifndef TRUE
#define TRUE (BOOL) 1
#endif
#ifndef FALSE
#define FALSE (BOOL) 0
#endif
#ifndef NULL
#define NULL (void *)(0)
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

/**********common structures*********************************************/
#define NDHEAPF 1200000L /* suggested length of store for heap (calfit)*/
#define MAXCAT 6  /* number of quanta in catalog file */
#define MAXQN 10  /* maximum number of quanta */
#define NDECDIP 6 /* number of digit pairs + 1 for idip */

typedef struct SXLINE_S {
  double cfrq;   /* calc. frequency   */
  double xfrq;   /* expt. frequency   */
  /*@dependent@*/ double *dnudp; /* freq. derivatives */
  float  xwt;    /* weight of line    */
  float  xerr;   /* error  of line    */
  int    linku;  /* link to next upper state */
  int    linkl;  /* link to next lower state */
  int    ibu;    /* upper state block number */
  int    ibl;    /* lower state block number */
  short  inu;    /* upper state index number */
  short  inl;    /* lower state index number */
  short  bln;    /* blend flag               */
  short  qn[2*MAXQN]; /* quantum number input     */
} SXLINE;
typedef struct SXLINE_S SXLINE;
typedef unsigned char bcd_t;
#define NEGBCD(ivbcd) (int)(ivbcd & 0x80)

#endif /* _CALPGM_TYPES_H_ */
