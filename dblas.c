#include <stdlib.h>
#include <math.h>
#include "cblas.h"
/*   Herbert M. Pickett, 22 March 1999 */
/*   Adapted from subset of BLAS routines in LINPAK */
CBLAS_INDEX cblas_idamax(const int n, const double *dx, const int incX)
{  /*  FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE. */
  double dmax, dtmp;
  int i, iret;
  if (n <= 1)
    return 0;
  dmax = fabs(dx[0]);
  iret = 0;
  for (i = 1; i < n; ++i) {
    dtmp = fabs(dx[i * incX]);
    if (dtmp > dmax) {
      dmax = dtmp;
      iret = i;
    }
  }
  return (CBLAS_INDEX) iret;
}                               /* idamax */

double cblas_dasum(const int N, const double *dx, const int incX)
{ /* TAKES THE SUM OF THE ABSOLUTE VALUES. */
  double dret = 0.0;
  int i;

  if (N <= 0)
    return 0.0;

  for (i = 0; i < N; ++i) {
    dret += fabs(dx[i * incX]);
  }

  return dret;
}    /* dasum */

void cblas_daxpy(const int N, const double da, const double *dx,
                 const int incX, double *dy, const int incY)
{  /* CONSTANT TIMES A VECTOR PLUS A VECTOR. */
  int i;

  if (N <= 0 || da == 0.0)
    return;

  for (i = 0; i < N; ++i) {
    dy[i * incY] += da * dx[i * incX];
  }
}  /* daxpy */

void cblas_dcopy(const int N, const double *dx, const int incX,
                 double *dy, const int incY)
{  /* COPIES A VECTOR, X, TO A VECTOR, Y. */
  int i;

  if (N <= 0)
    return;

  if (incX == 0) {
    /* Special case: copy the first element of dx to all elements of dy */
    double value = dx[0];
    for (i = 0; i < N; ++i) {
      dy[i * incY] = value;
    }
  } else {
    /* Normal case: copy each element */
    for (i = 0; i < N; ++i) {
      dy[i * incY] = dx[i * incX];
    }
  }
} /* dcopy */

double cblas_ddot(const int N, const double *dx, const int incX,
                  const double *dy, const int incY)
{   /*  FORMS THE DOT PRODUCT OF TWO VECTORS. */
  double dret = 0.0;
  int i;

  if (N <= 0)
    return 0.0;

  if (incX == incY && dx == dy) {
    /* Special case: dot product of a vector with itself */
    for (i = 0; i < N; ++i) {
      double val = dx[i * incX];
      dret += val * val;
    }
  } else {
    /* General case: dot product of two different vectors */
    for (i = 0; i < N; ++i) {
      dret += dx[i * incX] * dy[i * incY];
    }
  }

  return dret;
}                               /* ddot */

double cblas_dnrm2(const int n, const double *dx, const int incX)
{
  /* calculates the length of vector x */
  double dsum, xabs, xmax, r;
  int i, nsum;

  if (n <= 0)
    return 0.0;

  nsum = 0;
  dsum = 1.0;
  xmax = 0.0;

  for (i = 0; i < n; ++i) {
    xabs = fabs(dx[i * incX]);

    if (xabs > xmax) {
      if (nsum != 0) {
        r = xmax / xabs;
        dsum = 1.0 + (dsum * r) * r;
      }
      xmax = xabs;
      ++nsum;
    } else if (xabs != 0.0) {
      r = xabs / xmax;
      dsum += r * r;
      ++nsum;
    }
  }

  if (nsum > 1)
    xmax *= sqrt(dsum);

  return xmax;
}                               /* dnrm2 */

void cblas_drot(const int N, double *dx, const int incX,
                double *dy, const int incY, const double c, const double s)
{                               /* APPLIES A PLANE ROTATION. */
  double dtmpx, dtmpy;
  int i;

  if (N <= 0)
    return;

  for (i = 0; i < N; ++i) {
    dtmpx = dx[i * incX];
    dtmpy = dy[i * incY];
    dx[i * incX] = dtmpx * c + dtmpy * s;
    dy[i * incY] = dtmpy * c - dtmpx * s;
  }
}                               /* drot */

void cblas_drotg(double *sa, double *sb, double *c, double *s)
{
  double asa, asb, r, t;
  asa = (*sa);
  if (asa < 0.)
    asa = -asa;
  asb = (*sb);
  if (asb < 0.)
    asb = -asb;
  if (asa > asb) {              /* c > 0.707 */
    t = (*sb) / (*sa);
    r = sqrt(1. + t * t);
    *sa *= r;
    *c = 1. / r;
    *s = *sb = (*c) * t;
  } else if (asb != 0.) {       /* s > 0.707 */
    t = (*sa) / (*sb);
    r = sqrt(1. + t * t);
    *s = 1. / r;
    *c = (*s) * t;
    *sa = (*sb) * r;
    *sb = 1.;
    if ((*c) != 0.)
      *sb /= (*c);
  } else {
    *sa = *s = 0.;
    *sb = *c = 1.;
  }
}
void cblas_dscal(const int N, const double da, double *dx, const int incX)
{  /* SCALES A VECTOR BY A CONSTANT. */
  int i;

  if (N <= 0 || da == 1.0)
    return;

  for (i = 0; i < N; ++i) {
    dx[i * incX] *= da;
  }
}                               /* dscal */

void cblas_dswap(const int N, double *dx, const int incX,
                 double *dy, const int incY)
{  /* INTERCHANGES TWO VECTORS. */
  double dtmp;
  int i;

  if (N <= 0)
    return;

  for (i = 0; i < N; ++i) {
    dtmp = dx[i * incX];
    dx[i * incX] = dy[i * incY];
    dy[i * incY] = dtmp;
  }
}  /* dswap */
