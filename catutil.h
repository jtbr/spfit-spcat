#ifndef CATUTIL_H
#define CATUTIL_H

#ifdef __cplusplus
extern "C" {
#endif

int readqn(const char *qnstr, short *iqn, const int n);
void gupfmt(int igup, /*@out@*/ char *sgup);
int pcard(const char *card, /*@out@*/ double *val, const int nval,
          /*@null@*/ const int *fmtlen);

#ifdef __cplusplus
}
#endif

#endif /* CATUTIL_H */
