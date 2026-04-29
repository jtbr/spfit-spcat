/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include <algorithm>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CalCat.hpp"
#include "splib/calpgm_types.h"
#include "splib/blas_compat.h"
#include "splib/ulib.h"
#include "common/CalError.hpp"

int CalCat::qnfmt(short *iqu, int nqn, char *sqn)
{ /* subroutine to do quick conversion of quantum integers to characters */
  static int czero = (int)'0';
  int i, k, ix;
  char tqn1, tqn2;

  for (k = 0; k < nqn; ++k) {
    i = iqu[k];
    if (i < 0) {
      i = -i;
      ix = i / 10;
      i -= ix * 10;
      tqn1 = (char)(i + czero);
      if (ix == 0) {
        tqn2 = '-';
      } else if (ix < 27) {
        tqn2 = (char)(ix + ((int)'a' - 1));
      } else {
        tqn1 = tqn2 = '*';
      }
    } else {
      ix = i / 10;
      i -= ix * 10;
      tqn1 = (char)(i + czero);
      if (ix == 0) {
        tqn2 = ' ';
      } else if (ix < 10) {
        tqn2 = (char)(ix + czero);
      } else if (ix < 36) {
        tqn2 = (char)(ix + ((int)'A' - 10));
      } else {
        tqn1 = tqn2 = '*';
      }
    }
    ix = k + k;
    sqn[ix] = tqn2;
    sqn[ix + 1] = tqn1;
  }
  return 0;
}

int CalCat::simt(int isiz, int jsiz, double *s, double *t, double *u, double *wk)
{
  int i, j;
  double *sij, *si, *tc, *uc;
  /*    subroutine to do similarity transform */
  /*      S' = transpose(T) * S * U */
  if (isiz > 1) { /* do transform on the left */
    sij = s;
    for (j = 0; j < jsiz; ++j) {
      dcopy(isiz, sij, 1, wk, 1);
      tc = t;
      for (i = 0; i < isiz; ++i) {
        *sij = ddot(isiz, tc, 1, wk, 1);
        tc += isiz;
        ++sij;
      }
    }
  }
  if (jsiz > 1 && u != NULL) { /* do transform on the right */
    si = s;
    for (i = 0; i < isiz; ++i) {
      dcopy(jsiz, si, isiz, wk, 1);
      sij = si;
      uc = u;
      for (j = 0; j < jsiz; ++j) {
        *sij = ddot(jsiz, wk, 1, uc, 1);
        uc += jsiz;
        sij += isiz;
      }
      ++si;
    }
  }
  return 0;
}

void CalCatOutput::sort_cat_lines()
{
    // Mirrors sortn's dokey=FALSE comparison (bufcmp):
    //   primary key  — cols  0–12  (frequency, right-justified fixed-width)
    //   secondary key — cols 55–94 (quantum numbers)
    // Both compared lexicographically; for positive fixed-width fields this
    // equals numeric order.  stable_sort matches sortn's address-difference
    // tie-breaker (original order preserved for equal lines).
    std::stable_sort(cat_lines.begin(), cat_lines.end(),
                     [](const std::string &a, const std::string &b) {
                         int cmp = a.compare(0, 13, b, 0, 13);
                         if (cmp != 0) return cmp < 0;
                         // secondary key: cols 55 onwards
                         const size_t sec = 55;
                         bool a_has = a.size() > sec, b_has = b.size() > sec;
                         if (!a_has && !b_has) return false;
                         if (!a_has) return true;
                         if (!b_has) return false;
                         return a.compare(sec, std::string::npos,
                                          b, sec, std::string::npos) < 0;
                     });
}

