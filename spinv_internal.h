/* shared internal header file for spinv functionality */
#ifndef _SPINV_INTERNAL_H_
#define _SPINV_INTERNAL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdbool.h>
#include <ctype.h>
#ifndef BOOL
#define BOOL bool
#endif
#ifndef TRUE
#define TRUE (BOOL)1
#endif
#ifndef FALSE
#define FALSE (BOOL)0
#endif
#ifndef NULL
#define NULL (void *)(0)
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif
typedef unsigned char bcd_t;

#define ODD(i)  (((int)(i) & 1) != 0) /* Check if integer i is odd */
#define EVEN(i) (((int)(i) & 1) == 0) /* Check if integer i is even */
#define ODD2(i) (((int)(i) & 2) != 0) /* Check if bit 1 of integer i is set (i.e., i is 2 or 3 mod 4) */
#define TEST(i,mask) (((int)(i) & (mask)) != 0) /* Check if any bits in mask are set in integer i */
#define MOD(i,n) ((int)(i) % (int)(n)) /* Integer modulo operation */
#define C0 '\0' /* Null character */
#define HAM_DEBUG 0 /* Debug flag for Hamiltonian calculations, currently off */
#define MAXINT  12  /* Maximum number of dipole types for icase in setint */
#define MAXVIB  999 /* Maximum number of vibrational states */
#define MAXSPIN 9   /* Maximum number of spins */
#define MAXN_DIRCOS 359 /* MAX K in DIRCOS, also used as general upper K limit */
#define MAXVDEC 3   /* Maximum number of decimal digits for vibrational quanta in BCD */
#define NDECPAR 5   /* Number of digit pairs for idpar (parameter ID) not including vibrational part */
#define NSBCD   (2 * (NDECPAR + MAXVDEC + 1)) /* Size of BCD string buffer: (pairs + vib_pairs + sign_pair)*2chars_per_pair */
#define MAXII   20  /* Maximum value for 2*I (spin quantum number), so I <= 9.5 */
#define MAXNS   20  /* Maximum size of spin coupling arrays like lscom, iscom, jscom, (2 * MAXSPIN + 2) */

/* SPAR flags bits for parameter properties */
#define MNSQ     2  /* Bit for operator which is different from previous by only N*(N+1) */
#define MCOS_OK  4  /* Bit for operator which uses cosine (K-dependent part) from previous parameter */
#define MNOUNIT  8  /* Bit for operator which is not a unit matrix (K-dependent part is not constant) */
#define MMASK   (0x0fe0) /* Mask to clear bits 0-4 of SPAR flags (these are MNSQ, MCOS_OK, MNOUNIT, etc.) */
#define MIDEN   16  /* Bit for Identity operator under alpha symmetry (for I_tot) */

/* SPAR/DIRCOS loff/operator property flags */
#define MODD  0x040 /* Bit for operator that is inherently imaginary (e.g. Px) before phase conventions */
#define MSYM2 0x080 /* Bit for operator that is (I_alpha - I_-alpha) type for I_tot symmetry */
#define MLZ   0x100 /* Bit for Lz operator component */
#define MFCM  0x200 /* Bit for FC (Fourier Coefficient) being negative (sine series) */
#define MIMAG 0x400 /* Bit for operator that becomes effectively imaginary after phase conventions/symmetry */

/* SDIP flg bits for dipole properties */
#define MDIPI    8  /* Bit for imaginary electric dipole (after phase conventions) */
#define MINOQ    4  /* Bit for dipole operator which has NO Delta N = 0 (e.g., a commutator like [N^2, mu]) */
#define MELEC    2  /* Bit for electric dipole (vs magnetic) */

/* Indices for ixcom/jxcom arrays, which store quantum numbers for a block/state */
#define NDXCOM  8   /* Size of ixcom/jxcom arrays */
#define XDIM    0   /* 0: Dimension of sub-block */
#define XSYM    1   /* 1: Symmetry code for block (0=A, 1=Bx, 2=By, 3=Bz in D2) */
#define XNVAL   2   /* 2: 'N' rotational quantum number (excluding spin) */
#define XKBGN   3   /* 3: Beginning K quantum number for the Wang block */
#define XLVAL   4   /* 4: l (lambda or internal rotation) quantum number */
#define XVIB    5   /* 5: Vibrational state index */
#define XISYM   6   /* 6: Spin symmetry index (for I_tot basis) */
#define XIQN    7   /* 7: Additional spin symmetry quantum number information (for I_tot) */

/* Common Declarations */
typedef struct {  /* Vibrational state information */
  /*@dependent@*/ short *spt; /* Pointer to spin table (array of 2*spin values for each spin) */
  int nspstat;                /* Number of spin states (product of (2*I+1) for all spins) */
  short knmin[4];             /* Minimum K for each of the four D2 symmetries (A, Bx, By, Bz order may vary with oblate/prolate) */
  short knmax;                /* Maximum K for this vibrational state */
  short wt[5];                /* Statistical weights: wt[0-3] for D2 symmetries, wt[4] for special I_tot cases or axis indicator */
  short ewt[2];               /* E-symmetry statistical weights (e.g., for Cnv, Dnh) [0]=A-type E, [1]=B-type E for D6h etc. */
  short lvqn;                 /* l quantum number (or torsional state index, often |l|) */
  short lvupper;              /* Flags for l-doubling pair: bit 0 for upper state of pair, bit 1 if l-state is |K|=0 but effectively |l|=1, bits 2+ are 2*|l| or 2*torsional_index */
  short gsym;                 /* Symmetry group number (e.g. 2 for D2, 3 for C3v, etc.), odd if I_tot basis active for this state */
  short nqn;                  /* Number of spin quantum numbers for this state ( F1, F2...F), updated in setopt */
} SVIB;

typedef struct str_spar {  /* Local parameter information (parsed from input .par file) */
  /*@null@*/ /*@owned@*/ struct str_spar *next; /* Pointer to next parameter in a linked list (for a given vib. state and symmetry) */
  double zfac;          /* Overall scaling factor for the parameter, includes C6jj, C3jj products etc. */
  int ip;               /* Index of this parameter in the input 'par' array */
  unsigned int ipsym;   /* Packed symmetry and vibrational state identifier for fast lookup */
  short flags;          /* Bit flags describing operator properties (MNSQ, MCOS_OK, etc.) */
  signed char ksq;      /* Power of K^2 (N_z^2) */
  signed char fc;       /* Fourier coefficient index (for internal rotation or similar expansions) */
  signed char kavg;     /* Specific K_avg value if parameter is sampled at fixed K (NOFC=1), or K for 8x type operators */
  signed char msi1;     /* Index of the first spin involved (1-based, 0 for N) */
  signed char msi2;     /* Index of the second spin involved (1-based, 0 if none, -1 for N.N in commutator) */
  signed char mldel;    /* Change in l quantum number (Delta l) */
  unsigned char njq;    /* Power of N(N+1) operator */
  unsigned char mln;    /* Tensor order in N (rotational angular momentum) */
  unsigned char mld;    /* Tensor order of direction cosine operator part */
  unsigned char mkdel;  /* Change in K quantum number (Delta K) */
  unsigned char mins;   /* Power of N.S operator (S is the first spin) */
  unsigned char msznz;  /* Flag/type for S_z N_z or related operators */
  unsigned char euler;  /* Euler series flag/index (0 if not Euler, else type of Euler term) */
  unsigned char alpha;  /* Symmetry component for I_tot basis operators */
} SPAR;
typedef /*@null@*/ /*@owned@*/ SPAR *PSPAR; /* Pointer to SPAR struct */

typedef struct struct_ssp { /* Spin state pattern structure */
  /*@null@*/ /*@owned@*/ /*@reldef@*/ struct struct_ssp *next; /* Pointer to next SSP in a linked list */
  /*@notnull@*//*@owned@*/ short *sspt; /* Pointer to the actual spin table (array of 2*spin values, then coupling info) */
  int ssize;         /* Total number of spin states for this combination of spins (product of 2I+1) */
  int nitot;         /* I_tot symmetry type (e.g., 3 for C3v) for this spin combination, 0 if not I_tot */
} SSP;

typedef struct { /* Local dipole information (parsed from input .int file) */
  double fac;    /* Dipole moment value or scaling factor */
  short  flg;    /* Bit flags describing dipole operator properties (MELEC, MODD, etc.) */
  signed char kd;   /* Change in K quantum number (Delta K) */
  signed char ld;   /* Tensor order of direction cosine part of dipole operator */
  signed char ldel; /* Change in l quantum number (Delta l) */
  signed char fc;   /* Fourier coefficient index for dipole moment */
  signed char kavg; /* Specific K_avg value if dipole is sampled at fixed K (NOFC=1) */
} SDIP;

typedef struct
{                      /* Global settings and state variables */
  int mxspin;          /* Maximum value of 2*I for any spin encountered */
  int idiag;           /* Diagonalization/sorting option (0-5, or -1 for no diag) */
  int nvib;            /* Number of vibrational states defined */
  int nfit;            /* Number of parameters being fitted (positive IDPAR) */
  int nbkpj;           /* Number of unique (Vib,Sym,SpinPattern) blocks per F value */
  int ixz;             /* Bitmask for interaction restrictions (Delta N, Delta J, Delta F1, etc.) */
  int nqnn;            /* Number of basic rotational/rovib quantum numbers (e.g., N,Ka,Kc or N,K for sym top) */
  int nqn;             /* Total number of quantum numbers in output format (including vibrations and spins) */
  int maxqn;           /* Maximum number of QN fields to actually fill in output */
  int vibfac;          /* Factor for packing vibrational quanta (9, 99, or 999) */
  int parinit;         /* Flag/counter indicating if parameter structures (spar_head) have been initialized / number of vib states initialized */
  int maxblk;          /* Maximum number of sub-blocks (Wang blocks * spin states) in any F block */
  int nitot;           /* I_tot symmetry type (e.g. 3 for C3, 4 for C4 etc.), 0 if no I_tot active */
  int vibdec;          /* Number of decimal digits used for vibrational quanta in BCD parameter IDs */
  int esymdec;         /* Factor for E-symmetry weights (100 or 1000) */
  int msmask;          /* Mask to extract vibrational index from packed (mvs) integer */
  int nqn0;            /* Base number of quantum numbers before spins/vib added for format */
  int iqfmt0;          /* Base QN format code before spin/vib details */
  int maxwt;           /* Max number of distinct I_tot symmetry weight components (e.g. 2 for D6h E1/E2) */
  int stdphase;        /* Standard phase convention choice (bitmask for x,y,z components) */
  int phasemask;       /* Mask of symmetry components that have determined the standard phase */
  int g12;             /* G12 symmetry flag for Fourier coefficients (alternating sign with K) */
  unsigned int msshft; /* Bit shift to apply to spin state index before packing with vibrational index and symmetry */
  BOOL lsym;           /* True if l-doubling symmetry needs special handling for K=0 states (asymmetric top QN) */
  BOOL esym;           /* True if any E-symmetry states are present (triggers more complex weight handling) */
  BOOL oblate;         /* True if molecule is treated as oblate rotor */
  BOOL vibfmt;         /* True if vibrational quantum number is included in the output QN format */
  BOOL newlz;          /* True if Lz operator is defined explicitly in parameter IDs */
  BOOL nofc;           /* True if FC field in parameter ID specifies K_avg rather than Fourier coefficient index */
} GLOB;

typedef struct {   /* Cached data for getqn function to speed up repeated calls for the same block */
  int cblk;       /* Cached block number */
  int cnblk;      /* Cached number of sub-blocks in this block */
  int csblk;      /* Cached current sub-block index being processed */
  int cff;        /* Cached 2*F value for this block */
  int cwt[5];     /* Cached statistical weights for this block */
} GETQN;

/* Function Declarations */
struct SpinvContext;

/* Public interface functions, formerly in calpgm.h, now encapsulated in SpinvEngine */
int hamx(struct SpinvContext *ctx, const int iblk, const int nsiz, const int npar, const bcd_t *idpar,
         const double *par, /*@out@*/ double *egy, /*@out@*/ double *t,
         /*@out@*/ double *dedp, /*@out@*/ double *pmix, const BOOL ifdump);
int setint(struct SpinvContext *ctx, FILE *lu, /*@out@*/ BOOL *ifdiag, /*@out@*/ int *nsav, const int ndip,
           bcd_t *idip, /*@out@*/ int *isimag);
int intens(struct SpinvContext *ctx, const int iblk, const int isiz, const int jblk, const int jsiz,
           const int ndip, const bcd_t *idip, const double *dip,
           /*@out@*/ double *s);
int getqn(struct SpinvContext *ctx, const int iblk, const int indx, const int maxqn, /*@out@*/ short *iqn,
          /*@out@*/ int *idgn);
int setopt(struct SpinvContext *ctx, FILE *lu, /*@out@*/ int *nfmt, /*@out@*/ int *itd,
           /*@out@*/ int *nbcdpar, /*@out@*/ char *namfil);
int setfmt(struct SpinvContext *ctx, /*@out@*/ int *iqnfmt, int nfmt);
int setblk(struct SpinvContext *ctx, FILE *lu, const int npar, bcd_t *idpar, const double *par,
           int *nblkpf, int *negy);

/* internal functions */

int dclr(const int n1, const int n2, double *vec, const int ix);
int specop(const int neuler, BOOL *newblk, int *nsqj, int *ikq,
                  const int ksi, const int ksj, const int ni, const int nj,
                  const int ncos, double *wk, const short *ix,
                  const short *jx, const double par);
int specfc(struct SpinvContext* ctx, const int ifc, const int iv, const int jv,
                  const int kdel, const int ksi, const int ksj,
                  const int ncos, double *wk, const short *ix,
                  const short *jx);
int sznzfix(struct SpinvContext* ctx, const int sznz, const int ni, const int nj, int *ixcom,
                  int *jxcom, int *iscom, int *jscom);
int sznzop(struct SpinvContext* ctx, const int ni, const int nj, const int ksi, const int ksj,
                  const int *iscom, const int *jscom, const int ncos,
                  double *wk, const short *ix, const short *jx);
unsigned int blksym(const int *ixcom, const int *jxcom);
int ordham(const int nn, short *mask, double *egy, const short *isblk,
                  short *iswap);
int fixham(const int ndm, const int nn, double *t, double *egy,
                  double *p, const short *iswap);
BOOL kroll(const int nsizd, double *t, const int nsblk,
                  const short *sbkptr, const short *kmin);
int bestk(const int ndm, const int nsize, short *iqnsep, short *ibkptr,
                  short *itau, short *idx, double *t, double *egy, double *pmix,
                  double *wk);
int getqs(struct SpinvContext* ctx, const int mvs, const int iff, const int nsiz, const int kbgn,
                  int *ixcom, /*@out@*/ int *iscom, /*@out@*/ int *iv);
int idpars(SPAR * pspar, /*@out@*/ int *ksq, /*@out@*/ int *itp,
                  /*@out@*/ int *l, /*@out@*/ int *ld, /*@out@*/ int *kdel,
                  /*@out@*/ int *ins, /*@out@*/ int *si1, /*@out@*/ int *si2,
                  /*@out@*/ int *sznz, /*@out@*/ int *ifc, /*@out@*/ int *alpha,
                  /*@out@*/ int *ldel,/*@out@*/ int *kavg);
int getll(struct SpinvContext *ctx, const int llf, const int ld, const int ln, const int kd, /* signature was wrong in original code */
                  const int si1, const int si2, int *lscom, const int *iscom,
                  const int *jscom);
int getmask(struct SpinvContext* ctx, const int *xbra, const int *xket, const int kd, const int ldel,
                  const int loff, const int alpha);
double rmatrx(const int ld, const int lv, const int *ixcom,
                  const int *jxcom);
int symnsq(struct SpinvContext* ctx, const int inq, const int ins, const int *iscom,
                  const int *jscom, double *z);
int symksq(const int ikq, const int ksi, const int ksj, const int n,
                  double *wk, short *ix, short *jx);
int dpmake(const int nsize, double *dp, const double *t,
                  const int n, const double *wk, const short *ix,
                  const short *jx, const int isunit);
int pasort(struct SpinvContext* ctx, FILE * lu, const int npar, bcd_t *idpar,
                  const double *par);
int  idpari(struct SpinvContext* ctx, bcd_t *idval, int itp, /*@out@*/ SPAR * pspar);
int checksp(struct SpinvContext* ctx, const BOOL first, int si1, int si2, const short *iiv1,
                  const short *iiv2, double *zfac);
int  tensor(struct SpinvContext* ctx, double *z, const int *iscom, const int *jscom,
                  const int *lscom, const int *smap, int npair, int alpha);
int   setwt(SVIB * pvinfo, const int ivib, const int iax,
                  const int iwtpl, const int iwtmn, double vsym);
int   getwt(struct SpinvContext* ctx, SVIB * pvinfo, const int isym, const int iispin,
                  /*@out@*/ int *ivwt);
BOOL testwt(struct SpinvContext* ctx, SVIB *pvib1, SVIB *pvib2, int isym, int alpha);
int checkwt(struct SpinvContext* ctx, int *iwt, int *jwt);
int setgsym(struct SpinvContext* ctx, const int gsym);
int   getsp(struct SpinvContext* ctx, const bcd_t *ispnx, SVIB *pvinfo);
void  setsp(struct SpinvContext* ctx);
int   getqq(struct SpinvContext* ctx, const int iblk, /*@out@*/ int *f, /*@out@*/ int *iwtb,
                  /*@out@*/ short *sbkptr, /*@out@*/ short *kmin,
                  /*@out@*/ int *vs);
int dircos(struct SpinvContext* ctx, const int *xbra, const int *xket, const int ld,
                  const int kd, const int ncmax, /*@out@*/ double *direl,
                  /*@out@*/ short *ibra, /*@out@*/ short *iket,
                  const int ifup, const int loff, int mask,
                  /*@out@*/ int *isunit);
int ffcal(const int nff, const int kff, /*@out@*/ double *ff);

#ifdef __cplusplus
} // end extern "C"
#endif

#endif /* _SPINV_INTERNAL_H_ */