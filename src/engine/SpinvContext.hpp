#ifndef SPINV_CONTEXT_HPP
#define SPINV_CONTEXT_HPP

#ifdef __cplusplus
extern "C" {
#endif

#include "splib/calpgm_types.h"
#include "splib/blas_compat.h"
#include "splib/ulib.h"
#include "splib/cnjj.h"
#include "spinit.h"
#include "spinv_internal.h"

/* Encapsulation of previous global context */
struct SpinvContext {
  SVIB vinfo1; /* Default/first vibrational state info */
  SVIB *vinfo; /* Pointer to array of vibrational state info */
  PSPAR spar_head[MAXVIB]; /* Array of linked list heads for parameters, one per vibrational state (iv2, the lower state in an interaction) */
  short sptzero[2]; /* Default spin table for a molecule with no spins: nset=1 (N only), 2*S=0 */
  SSP ssp_head; /* Head of linked list for unique spin patterns */
  GETQN *cgetq;                  /* Two cache entries, presumably for upper and lower states of a transition */
  GETQN cgetv[2];
  SDIP dipinfo0;                 /* Default/first dipole info structure */
  SDIP *dipinfo;                 /* Pointer to array of dipole info structures */

  double zwk;                    /* Static zero for work array fallback */
  double spfac[MAXII];           /* Spin factors sqrt(I(I+1)(2I+1)) for quadrupole etc. */
  double spfac2[MAXII];          /* Spin factors for quadrupole, involving 1/sqrt((I-1/2)(I+3/2)) etc. */
  int zmoldv;                    /* Default/zero for moldv fallback */
  int zblkptr;                   /* Default/zero for blkptr fallback */
  int zivs;                      /* Default/zero for ivs fallback */
  int zipder;                    /* Default/zero for ipder fallback */
  int revsym[4];  /* Symmetry reversal map (e.g., Bx <-> Bz for oblate) */
  int isoddk[4];  /* Indicates if K is effectively odd for symmetry types A, Bx, By, Bz for direction cosines */
  int ixphase[4]; /* Scratch for phase choices, modified by glob.stdphase */
  int ipwr2[4];   /* Powers of 2 (1,2,4) used for symmetry component checking (bit masks) */
  int is_esym[MAXITOT];          /* Array indicating E-symmetry for I_tot components: 0=A/B, 1=E_a, -1=E_b */

  /* Angular momentum coupling arrays for a single operator */
  int lscom[MAXNS];  /* Tensor orders L for successive couplings in F = N+S+I1+I2... */
  int iscom[MAXNS];  /* 2*Angular momentum values for bra state (2N, 2J, 2F1...) */
  int jscom[MAXNS];  /* 2*Angular momentum values for ket state */
  int ismap[MAXNS];  /* Maps spin index to its position in iscom/jscom for tensor calculations */

  /* Common blocks for rotational/rovibrational quanta storage */
  int ixcom[NDXCOM]; /* Stores quantum numbers for the 'bra' side of a matrix element */
  int jxcom[NDXCOM]; /* Stores quantum numbers for the 'ket' side of a matrix element */

  /* Various global integer flags and counters */
  int itptr;  /* Index of I_tot in the spin coupling scheme */
  int itsym;  /* Index of the first spin summed into I_tot */
  int nspin;  /* Maximum number of spins encountered across all vibrational states */
  int nsqmax; /* Maximum power of N(N+1) encountered for any parameter */
  int ndmx;   /* Maximum Hamiltonian matrix dimension encountered / work array size */
  int ndmax;  /* Maximum Delta N for any interaction */
  int nddip;  /* Number of dipole parameters currently allocated in dipinfo */

  /* Static short variables, often default/zero values or fallbacks for unallocated pointers */
  short zidx;    /* Default/zero for idx fallback */
  short zjdx;    /* Default/zero for jdx fallback */
  short ziqnsep; /* Default/zero for iqnsep fallback */
  short zibkptr; /* Default/zero for ibkptr fallback */
  short zikmin;  /* Default/zero for ikmin fallback */

  GLOB glob;     /* Global variables structure */

  char sbcd[NSBCD]; /* Buffer for BCD (Binary Coded Decimal) to string conversion */

  /* Pointers to dynamically allocated arrays used throughout calculations */
  int *moldv;     /* Array storing packed (Vib,Sym,SpinPattern) identifiers for each block */
  int *blkptr;    /* Array of pointers to the start of each F-block group in moldv */
  int *ipder;     /* Array mapping parameter index to its derivative column, or -ve for constrained */
  int *ivs;       /* Array storing packed (Vib,Sym,SpinPattern) for each sub-block */
  double *wk;     /* Main work array for Hamiltonian elements, direction cosines, etc. */
  short *idx;     /* Row indices for sparse matrix elements (Hamiltonian or dipole) */
  short *jdx;     /* Column indices for sparse matrix elements */
  short *iqnsep;  /* Array storing separation information for sorting/diagonalization, or K values for projection sort */
  short *ibkptr;  /* Array of pointers to the start of each sub-block (Wang block * spin) */
  short *ikmin;   /* Array of minimum K values for each sub-block */
};

#ifdef __cplusplus
} // extern "C"
#endif

#endif // SPINV_CONTEXT_HPP