#ifndef _BLAS_COMPAT_H_
#define _BLAS_COMPAT_H_

/* Short-name aliases for CBLAS functions used throughout the codebase.
   cblas.h is the vendored upstream Netlib CBLAS header. */
#ifdef __cplusplus
extern "C" {
#endif
#include "cblas.h"
#ifdef __cplusplus
}
#endif
#define idamax cblas_idamax
#define dasum  cblas_dasum
#define daxpy  cblas_daxpy
#define dcopy  cblas_dcopy
#define ddot   cblas_ddot
#define dnrm2  cblas_dnrm2
#define drot   cblas_drot
#define drotg  cblas_drotg
#define dscal  cblas_dscal
#define dswap  cblas_dswap

#endif /* _BLAS_COMPAT_H_ */
