#ifndef CATUTIL_H
#define CATUTIL_H

#ifdef __cplusplus
extern "C" {
#endif

int readqn(const char *qnstr, short *iqn, const int n);
void gupfmt(int igup, /*@out@*/ char *sgup);

#ifdef __cplusplus
}
#endif

#endif /* CATUTIL_H */
