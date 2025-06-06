/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in c, 22 March 1999 */
/*   24 March  1999: Fixed sin scaling in dircos */
/*   25 March  1999: Fixed diag = 2 error in call to ordham*/
/*   27 March  1999: Changed diag = 3 to take care of l-doubled states */
/*   30 March  1999: Fixed sin scaling bug */
/*   21 April  1999: Fixed bug at end of hamx for KROLL case */
/*   22 April  1999: Fixed bug in Itot option */
/*   25 May    1999: fix bug when largest number of spins is not default */
/*   23 July   1999: fix bug in SETOPT that returned incorrect IQNFMT */
/*    7 Sept.  1999: fix bug in GETQQ for that filled ivs with short */
/*   29 Dec.   1999: improve ORDHAM */
/*   12 April  2000: negate Hamiltonian for kappa > 0 on first pass */
/*   12 April  2000: separate ORDHAM and FIXHAM */
/*   12 April  2000: subtract constant from energy for more precision*/
/*    1 August 2000: special code for unit matrix */
/*   12 Sept.  2000: fix up phase for cross-multiplet interaction */
/*    9 Oct.   2000: change dipole type 10 to double commutator with N^2 */
/*   29 Nov.   2000: redefine idpars and spar to put more calc. in initial */
/*    5 Feb.   2001: define augmented Euler and Fourier series */
/*    7 July   2001: fix bug for l doublet pairs in setopt */
/*   26 Sept.  2001: fix bug in idpari for oblate 6x and 2x interactions */
/*   12 March  2002: improve selections for off-diagonal l in dircos */
/*   27 June   2002: increase vibration limit to 359 */
/*   29 July   2002: fix bug in fixham (DIAG = 2,3)*/
/*   18 Aug.   2003: code cleanup, @comment@ is for splint */
/*   24 Nov.   2003: implement code for higher symmetry IAX options */
/*    3 Feb    2004: implement code for higher symmetry spins */
/*   19 May    2004: implement code for bcd storage of idpar and idip */
/*   20 July   2004: fix bug in code for bcd storage of idip */
/*   19 Jan.   2005: extended Euler operator */
/*   28 May    2005: fix bug in getmask for operator between l-doubled state */
/*   27 Sept.  2005: fix bug for 6-fold symmetry B-weights in getqn & setopt */
/*   17 Oct.   2005: fix bug in getmask for operator within l-doubled state */
/*   12 Dec.   2005: fix bug Itot basis and include mixed basis */
/*   10 Sept.  2006: add other phase conventions */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calpgm.h"
#include "spinit.h"
#include "spinv_internal.h"

SVIB vinfo1; /* Default/first vibrational state info */
/*@owned@*/ SVIB *vinfo = &vinfo1; /* Pointer to array of vibrational state info */


PSPAR spar_head[MAXVIB]; /* Array of linked list heads for parameters, one per vibrational state (iv2, the lower state in an interaction) */

short sptzero[] = { 1, 0 }; /* Default spin table for a molecule with no spins: nset=1 (N only), 2*S=0 */

SSP ssp_head = { NULL, sptzero, 1 , 0}; /* Head of linked list for unique spin patterns */

struct {   /* Cached data for getqn function to speed up repeated calls for the same block */
  int cblk;       /* Cached block number */
  int cnblk;      /* Cached number of sub-blocks in this block */
  int csblk;      /* Cached current sub-block index being processed */
  int cff;        /* Cached 2*F value for this block */
  int cwt[5];     /* Cached statistical weights for this block */
} *cgetq, cgetv[2]; /* Two cache entries, presumably for upper and lower states of a transition */


SDIP dipinfo0; /* Default/first dipole info structure */
/*@owned@*/ SDIP *dipinfo = &dipinfo0; /* Pointer to array of dipole info structures */

double zero = 0.;    /* Static zero value */
double zwk = 0.;     /* Static zero for work array fallback */
double spfac[MAXII]; /* Spin factors sqrt(I(I+1)(2I+1)) for quadrupole etc. */
double spfac2[MAXII];/* Spin factors for quadrupole, involving 1/sqrt((I-1/2)(I+3/2)) etc. */
int zmoldv = 0;      /* Default/zero for moldv fallback */
int zblkptr = 0;     /* Default/zero for blkptr fallback */
int zivs = 0;        /* Default/zero for ivs fallback */
int zipder = 0;      /* Default/zero for ipder fallback */
int revsym[] = { 0, 3, 2, 1}; /* Symmetry reversal map (e.g., Bx <-> Bz for oblate) */
int isoddk[] = { 0, 1, 1, 0}; /* Indicates if K is effectively odd for symmetry types A, Bx, By, Bz for direction cosines */
int ixphase[] = {0, 1, 2, 3}; /* Scratch for phase choices, modified by glob.stdphase */
int ipwr2[]={0,1,2,4};        /* Powers of 2 (1,2,4) used for symmetry component checking (bit masks) */
int is_esym[MAXITOT];         /* Array indicating E-symmetry for I_tot components: 0=A/B, 1=E_a, -1=E_b */
/* Angular momentum coupling arrays for a single operator */
int lscom[MAXNS];  /* Tensor orders L for successive couplings in F = N+S+I1+I2... */
int iscom[MAXNS];  /* 2*Angular momentum values for bra state (2N, 2J, 2F1...) */
int jscom[MAXNS];  /* 2*Angular momentum values for ket state */
int ismap[MAXNS];  /* Maps spin index to its position in iscom/jscom for tensor calculations */
/* Common blocks for rotational/rovibrational quanta storage */
int ixcom[NDXCOM]; /* Stores quantum numbers for the 'bra' side of a matrix element */
int jxcom[NDXCOM]; /* Stores quantum numbers for the 'ket' side of a matrix element */
/* Various global integer flags and counters */
int itptr;         /* Index of I_tot in the spin coupling scheme */
int itsym;         /* Index of the first spin summed into I_tot */
int nspin;         /* Maximum number of spins encountered across all vibrational states */
int nsqmax;        /* Maximum power of N(N+1) encountered for any parameter */
int ndmx;          /* Maximum Hamiltonian matrix dimension encountered / work array size */
int ndmax;         /* Maximum Delta N for any interaction */
int nddip;         /* Number of dipole parameters currently allocated in dipinfo */
/* Static short variables, often default/zero values or fallbacks for unallocated pointers */
short szero = 0;   /* Static zero short value */
short zidx = 0;    /* Default/zero for idx fallback */
short zjdx = 0;    /* Default/zero for jdx fallback */
short ziqnsep = 0; /* Default/zero for iqnsep fallback */
short zibkptr = 0; /* Default/zero for ibkptr fallback */
short zikmin = 0;  /* Default/zero for ikmin fallback */

GLOB glob; /* Global variables structure */

char sbcd[NSBCD]; /* Buffer for BCD (Binary Coded Decimal) to string conversion */

/* Pointers to dynamically allocated arrays used throughout calculations */
/*@owned@*/ int *moldv;   /* Array storing packed (Vib,Sym,SpinPattern) identifiers for each block */
/*@owned@*/ int *blkptr;  /* Array of pointers to the start of each F-block group in moldv */
/*@owned@*/ int *ipder;   /* Array mapping parameter index to its derivative column, or -ve for constrained */
/*@owned@*/ int *ivs;     /* Array storing packed (Vib,Sym,SpinPattern) for each sub-block */
/*@owned@*/ double *wk;   /* Main work array for Hamiltonian elements, direction cosines, etc. */
/*@owned@*/ short *idx;   /* Row indices for sparse matrix elements (Hamiltonian or dipole) */
/*@owned@*/ short *jdx;   /* Column indices for sparse matrix elements */
/*@owned@*/ short *iqnsep;/* Array storing separation information for sorting/diagonalization, or K values for projection sort */
/*@owned@*/ short *ibkptr;/* Array of pointers to the start of each sub-block (Wang block * spin) */
/*@owned@*/ short *ikmin; /* Array of minimum K values for each sub-block */


/*****************************************************************************/
/**
 * @brief Calculates Hamiltonian matrix elements, energies, and derivatives for a specific block.
 * @param iblk Block number to calculate.
 * @param nsize Dimension of the Hamiltonian block.
 * @param npar Total number of parameters.
 * @param idpar Array of parameter identifiers (BCD format), idpar[0] is length.
 * @param par Array of parameter values.
 * @param egy Output array for eigenvalues (energies).
 * @param t Input/Output array for Hamiltonian matrix elements, becomes eigenvector matrix on output.
 * @param dedp Output array for derivatives of energies with respect to parameters (dE/dP).
 * @param pmix Output array for mixing coefficients or state labels (e.g. K for projection sort).
 * @param ifdump If TRUE, dump Hamiltonian matrix `t` without diagonalizing and store diagonal elements in `egy`.
 * @return int Status: 1 for success, -1 if idpar[0] is 0 (potentially indicating an issue or no parameters).
 */
int hamx(iblk, nsize, npar, idpar, par, egy, t, dedp, pmix, ifdump)
const int iblk, nsize, npar;
const BOOL ifdump;
const bcd_t *idpar;
const double *par;
double *egy, *t, *dedp, *pmix;
{
  /* .. PACKAGE FOR 98 INTERACTING VIBRATIONAL STATES WITH MULTI-SPIN */

  /*     calculate energies and derivatives for block IBLK */

  static double sqj[10]; /* Array to store powers of N(N+1) scaled by ZPAR */
  static int nold = -1;  /* Previous N value, to optimize N(N+1) calculation */
  SPAR *spar_now;        /* Current parameter structure being processed */
  short *itau, *vbkptr;  /* itau: array for tau quantum numbers for sorting; vbkptr: pointers for vibrational sub-blocks */
  double *pt, *dbar, pbar, sqnn, zpar, dn, ele, egy0; /* pt: temp pointer into matrix t; dbar: derivative accumulator for current operator; pbar: composite parameter value for current operator; sqnn: N(N+1); zpar: K-independent part of operator; dn: N as double; ele: element contribution; egy0: energy offset for diagonal operators */
#if HAM_DEBUG
  double *pscr; /* Scratch space for debugging Hamiltonian */
#endif
  long ndm;             /* Leading dimension of matrix t (nsize) */
  int ipar, ncos, ispn, jspn, nsqj, iiwt[5], ndmd, alpha, mkd; /* ipar: parameter index; ncos: number of non-zero direction cosine elements; ispn,jspn: packed spin state for bra/ket; nsqj: power of N(N+1); iiwt: cached weights for current block; ndmd: nsize+1 for diagonal indexing; alpha: Itot symmetry component; mkd: mask for dircos */
  int i, ii, lt, n_sub_block_size, ibase, jbase, kbgni, kbgnj, nsblk, idflags, ipbase; /* i,ii: loop counters; lt: tensor L order; n_sub_block_size: current sub-block size; ibase,jbase: base index for bra/ket sub-block; kbgni,kbgnj: K_begin for bra/ket; nsblk: number of sub-blocks in current F block; idflags: flags from SPAR struct; ipbase: base parameter index for constrained params */
  int kd, ld, nd, ni, nj, ivbase, jvbase, ivmin, ifc, iz, jz, npair, ldel; /* kd: Delta K; ld: dir.cos. L order; nd: N_bra - N_ket; ni,nj: N_bra, N_ket; ivbase,jvbase: vib. index for bra/ket; ivmin: min(ivbase,jvbase); ifc: Fourier coeff index; iz,jz: matrix indices; npair: number of spin couplings; ldel: Delta l */
  int si1, si2, iff, ijd, ikq, ins, neuler, sznz, ixx, jxx, kl, isgn, isunit; /* si1,si2: spin indices; iff: 2*F; ijd: ixx-jxx; ikq: K^2 power; ins: N.S power; neuler: Euler type; sznz: SzNz type; ixx,jxx: sub-block loop indices; kl: loff flags for dircos; isgn: sign factor for matrix element; isunit: flag if operator is unit matrix type */
  int kavg, nqnsep; /* kavg: K_average for specific K parameters; nqnsep: result from hdiag, number of separated states or error */
  unsigned int ivsym, lastvv, ivcmp; /* Packed vibrational and symmetry identifiers for fast parameter lookup */
  BOOL oldpar, isegy0, roll, parskp, isneg, first; /* oldpar: if previous param shared cosine part; isegy0: if current operator contributes to E0; roll: if K-roll detected; parskp: skip current parameter; isneg: if oblate hamiltonian needs negation; first: first pass (calculating H) vs second pass (derivatives) */
  BOOL newblk, firstpar; /* newblk: if new vibrational pair, re-eval Euler denoms; firstpar: first parameter for a given vib pair and symmetry */

  if (ndmx <= 0) { /* Check if work arrays are allocated */
    puts("working vectors not allocated");
    exit(EXIT_FAILURE);
  }
  dbar = &wk[ndmx]; /* dbar points to a section of the work array wk, used for accumulating derivatives for a single operator */
  ndm = nsize;      /* Store nsize in ndm (leading dimension for matrix t) */
  ndmd = nsize + 1; /* Used for indexing diagonal elements t[i+i*ndm] as t[i*ndmd] */
  isneg = (glob.oblate && !ifdump); /* Hamiltonian needs to be negated for oblate rotors if not just dumping */
  cgetv[0].cblk = 0; /* Invalidate cache for getqn for state 0 */
  cgetv[1].cblk = 0; /* Invalidate cache for getqn for state 1 */
  /*     get F and sub-block structure */
  itau = &ikmin[glob.maxblk]; /* itau points to an area after ikmin, used to store sorting keys (tau or K values) */
  nsblk = getqq(iblk, &iff, iiwt, ibkptr, ikmin, ivs); /* Get sub-block structure for current F block (iblk) */
  /* zero hamiltonian matrix t */
  dclr(nsize, nsize, t, 1);
  /* zero derivative matrix dedp */
  n_sub_block_size = glob.nfit; /* Number of parameters to fit */
  dclr(nsize, n_sub_block_size, dedp, 1);
  /* set up hamiltonian */
  oldpar = parskp = roll = firstpar = FALSE; /* Initialize flags */
  ncos = sznz = 0; isgn = 1; isunit = 0; kbgni = 0; /* Initialize operator properties */
  pbar = zpar = 0.; /* Initialize composite parameter and K-independent factor */
  first = TRUE;     /* Start with the first pass (Hamiltonian construction) */
  do { /* This outer loop runs twice: once for H (first=TRUE), once for dH/dP (first=FALSE) */
    /* loop on wang blocks (actually, pairs of sub-blocks for interactions) */
    lastvv = 0; newblk = TRUE; egy0 = 0.; /* lastvv tracks vib pair to re-init Euler denoms if needed */
    for (ixx = 0; ixx < nsblk; ++ixx) { /* Loop over 'bra' sub-blocks (Wang blocks * spin states) */
      ibase = ibkptr[ixx]; /* Starting row/col index for this 'bra' sub-block */
      n_sub_block_size = ibkptr[ixx + 1] - ibase; /* Size of this 'bra' sub-block */
      kbgni = ikmin[ixx]; /* Starting K for this 'bra' sub-block */
      ispn = ivs[ixx];    /* Packed (Vib,Sym,SpinPattern) for 'bra' */
      getqs(ispn, iff, n_sub_block_size, kbgni, ixcom, iscom, &ivbase); /* Get quantum numbers for 'bra' state */
      ni = ixcom[XNVAL]; /* N quantum number for 'bra' */
      /* set up itau for sorting (tau = Ka-Kc like proxy, or K for projection) */
      kd = (ni + ixcom[XSYM] + 1) & 1; /* Parity component for tau */
      kd += kbgni + kbgni;             /* K component for tau */
      kd = kd + kd;                    /* Scale it up */
      if (ivbase != ixcom[XVIB]) /* Check if it's the upper state of an l-doublet pair based on how getqs returns XVIB relative to the actual state index ivbase */
        ++kd; /* Adjust tau for upper l-doubled state */
      itau[ixx] = (short) kd; /* Store sorting key for this sub-block */

      if (ni != nold) {  /* If N changed, re-calculate powers of N(N+1) */
        nold = ni;
        dn = (double) ni;
        sqnn = dn * (ni + 1); /* N(N+1) */
        sqj[0] = 1.0; /* Zeroth power */
        for (i = 1; i <= nsqmax; ++i) { /* Calculate higher powers */
          sqj[i] = sqj[i - 1] * sqnn;
        }
      }
      for (jxx = 0; jxx <= ixx; ++jxx) { /* Loop over 'ket' sub-blocks (only jxx <= ixx for lower triangle + diagonal) */
        ijd = ixx - jxx; /* Difference in sub-block indices */
        jbase = ibkptr[jxx]; /* Starting row/col index for this 'ket' sub-block */
        n_sub_block_size = ibkptr[jxx + 1] - jbase; /* Size of this 'ket' sub-block */
        kbgnj = ikmin[jxx]; /* Starting K for this 'ket' sub-block */
        jspn = ivs[jxx];    /* Packed (Vib,Sym,SpinPattern) for 'ket' */
        getqs(jspn, iff, n_sub_block_size, kbgnj, jxcom, jscom, &jvbase); /* Get quantum numbers for 'ket' state */
        nj = jxcom[XNVAL]; /* N quantum number for 'ket' */
        nd = ni - nj;      /* Delta N */
        if (nd > ndmax && ndmax >=0) continue; /* Skip if Delta N too large (ndmax < 0 allows all Delta N) */

        ivmin = (ivbase < jvbase) ? ivbase : jvbase; /* Minimum of the two vibrational indices */
        /* Create a packed identifier for the vibrational pair and overall block symmetry */
        ivsym = (unsigned int) (ivbase + jvbase + ivmin * glob.vibfac) << 2;
        ivcmp = ivsym + 3; /* Upper bound for checking Euler parameters which might be only vv' specific */
        ivsym += blksym(ixcom, jxcom); /* Add block symmetry (A,Bx,By,Bz) */

        if (lastvv != ivcmp) { /* If vibrational pair changed, Euler denominators might need re-init */
          newblk = TRUE; lastvv = ivcmp;
        }
        isegy0 = FALSE; firstpar = TRUE; /* isegy0: current op contributes to E0; firstpar: first param for this (v,v',sym) */

        for (spar_now = spar_head[ivmin]; TRUE; spar_now = spar_now->next) { /* Loop over parameters relevant to this vibrational pair */
          idflags = 0;
          if (spar_now != NULL) { /* Check if current parameter matches current ivsym/ivcmp */
            ivcmp = spar_now->ipsym;
            if (ODD(spar_now->euler)) { /* Euler denominator parameter (applies to a range of symmetries) */
              if (ivcmp > lastvv) continue; /* Parameter is for a later vv' pair */
              if (ivcmp == lastvv)
                idflags = (int) spar_now->flags; /* Matched Euler denominator parameter */
            } else { /* Normal operator parameter */
              if (ivcmp > ivsym) continue; /* Parameter is for a later symmetry or vv' */
              if (ivcmp == ivsym)
                idflags = (int) spar_now->flags; /* Matched operator parameter */
            }
          }
          /* This block processes the previously accumulated operator (if any) */
          if (!TEST(idflags, MCOS_OK)) { /* If current param doesn't share K-dep part, or end of list */
            if (oldpar) { /* If there was a previous operator sharing K-dep parts */
              /*  add diagonal or sub-diagonal to Hamiltonian */
              /*  using composite parameter value PBAR */
              if (isegy0) { /* If it's an E0 type term */
                egy0 += pbar;
              } else { /* Normal Hamiltonian term */
                if (isneg) /* Negate for oblate rotor representation */
                  pbar = -pbar;
                if (isunit != 0) { /* If operator is unit matrix type (K-independent part is 1) */
                  iz = idx[0]; /* Only one set of indices for unit matrix */
                  jz = jdx[0];
                  pt = &t[jz * ndm]; /* Pointer to start of column jz */
                  for (i = 0; i < ncos; ++i) { /* ncos is actual # elements for unit matrix */
                    pt[iz] += pbar; /* Add to H[iz, jz] (actually H[row_idx[i], col_idx[i]])*/
                    iz += ndmd; /* Move to next element assuming it's effectively on a diagonal of a submatrix if isunit represents block diagonal structure */
                                /* This indexing for unit matrix seems specific; idx[0], jdx[0] are start, ncos is length of this unit block along diagonal */
                  }
                } else { /* Operator has K-dependent part stored in wk array */
                  for (i = 0; i < ncos; ++i) {
                    iz = idx[i]; /* Row index for this element */
                    jz = jdx[i]; /* Column index for this element */
                    t[iz + jz * ndm] += pbar * wk[i]; /* Add PBAR * K_dependent_factor */
                  }
                }
              }
              oldpar = FALSE; /* Reset flag */
            }
            if (idflags == 0) break; /* End loop over parameters for this (v,v',sym) if no more matching params */
            ncos = -1; /* Signal that K-dependent part (wk, idx, jdx) needs recalculation */
          } else { /* Current param shares K-dep part with previous (TEST(idflags, MCOS_OK) is true) */
            if (ncos == 0) continue; /* If previous K-dep part was zero, skip this param too */
          }
          if (spar_now == NULL) break; /* Should have been caught by idflags == 0 if spar_now became NULL */

          ipar = spar_now->ip; /* Index of the parameter value in main 'par' array */
          nsqj = (int) spar_now->njq; /* Power of N(N+1) */

          if (TEST(idflags, MNSQ) && nd == 0) { /* If param differs from prev only by N(N+1) and DeltaN=0 */
            if (parskp) continue; /* Skip if the 'mother' operator (without this N(N+1) term) was zero */
          } else { /* New type of operator, or one not differing by just N(N+1) */
            parskp = TRUE; /* Assume initially this new operator type might be zero */
            if (sznz < 0) /* If SzNz type needs finalization (from previous param processing) */
              sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom); /* Reset N values in ixcom/jxcom if changed by SzNz */
            /* Parse the BCD parameter ID into its components */
            kl = idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins,
                        &si1, &si2, &sznz, &ifc, &alpha, &ldel, &kavg);
            if (sznz > 0) { /* If this is an SzNz type operator */
              sznz = sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom); /* Modify N in ixcom/jxcom, returns new sznz state */
              if (sznz == 0) continue;    /* Operator invalid for these N values, get next parameter */
            }
            if (kavg > 0 && kavg > ni && kavg > nj) continue; /* If specific K_avg is outside range of N_bra, N_ket */

            mkd = getmask(ixcom, jxcom, kd, ldel, kl, alpha); /* Get selection mask for dircos based on symmetry */
            if (mkd == 0) continue; /* Operator forbidden by symmetry/selection */

            npair = getll(0, ld, lt, kd, si1, si2, lscom, iscom, jscom); /* Get spin coupling tensor orders */
            if (npair < 0) continue; /* Tensor coupling rules violated */

            if (ncos < 0) { /* If K-dependent part needs recalculation */
              isunit = kavg; /* kavg might indicate a unit matrix type if > 0 and nofc=1 */
              /* Calculate direction cosine matrix elements (K-dependent part) */
              ncos = dircos(ixcom, jxcom, ld, kd, ndmx, wk, idx, jdx,
                            ijd, kl, mkd, &isunit);
              if (ncos <= 0) {
                if (ncos == 0) continue; /* No non-zero elements, get next parameter */
                printf("DIRCOS WORKING VECTOR TOO SHORT IN HAMX BY %d\n", -ncos); /* Error if ncos is negative */
                exit(EXIT_FAILURE);
              }
              if (TEST(idflags, MNOUNIT)) /* If flag says it's not a unit matrix, ensure isunit is 0 */
                isunit = 0;
              isgn = 1;
              if (isunit != 0 && wk[0] < 0.) /* For unit matrix type, if scaling factor is negative */
                isgn = -1;

              if (neuler != 0) { /* If Euler series operator */
                /* Apply Euler series transformation to wk, idx, jdx */
                /* par[ipar] is the Euler denominator 'a' or 'b' for this specific neuler type */
                ncos = specop(neuler, &newblk, &nsqj, &ikq, kbgni, kbgnj,
                              ni, nj, ncos, wk, idx, jdx, par[ipar]);
                if (ncos == 0) continue; /* Euler op became zero */
              }
              if (ifc != 0) /* If Fourier coefficient operator */
                specfc(ifc, ivbase, jvbase, kd, kbgni, kbgnj, ncos,
                       wk, idx, jdx); /* Apply Fourier scaling to wk */
              if (sznz != 0) /* If SzNz type operator */
                sznzop(ni, nj, kbgni, kbgnj, iscom, jscom, ncos, wk, idx, jdx); /* Apply SzNz scaling to wk */
              if (ikq > 0) { /* If K^2 dependent operator */
                ncos = symksq(ikq, kbgni, kbgnj, ncos, wk, idx, jdx); /* Apply K^2 scaling */
                if (ncos == 0) continue;
              }
              /* Adjust indices if not in the first sub-block of the matrix */
              if (ibase != 0) {
                for (i = 0; i < ncos; ++i) idx[i] = (short) (idx[i] + ibase);
              }
              if (jbase != 0) {
                for (i = 0; i < ncos; ++i) jdx[i] = (short) (jdx[i] + jbase);
              }

              if (first) { /* First pass (Hamiltonian construction) */
                isegy0 = FALSE;
                pbar = 0.; /* Initialize composite parameter for this K-dependent operator */
                if (isunit != 0 && ixx == jxx) { /* If diagonal unit operator */
                  oldpar = TRUE; /* This K-dependent part can be reused */
                  if (ixx == 0) { /* If it's in the very first Wang block, it's an E0 term */
                    isegy0 = TRUE;
                  } else if (firstpar){ /* If it's the first param for this vib/sym, subtract previous E0 */
                    pbar = -egy0;
                    firstpar = FALSE;
                  }
                }
                if (ipar < 0) continue; /* Should not happen if idflags was set from spar_now */
              } else { /* Second pass (derivatives dH/dP) */
                if (ipar < 0) continue;
                /* Calculate <psi| (dH_op / dP_current) |psi> using current eigenvectors t */
                /* dbar stores this for each state, effectively (d H_matrix_element / dP) transformed by eigenvectors */
                dpmake(nsize, dbar, t, ncos, wk, idx, jdx, isunit);
              }
            } /* End K-dependent part calculation (if ncos < 0) */

            if (si2 < 0 && nd == 0 && ni ==0) continue;  /* Avoids issues with commutator with N*N when N=0 */

            /* Calculate K-independent part of operator matrix element */
            zpar = rmatrx(ld, lt, ixcom, jxcom); /* N-dependent corrections for reduced matrix elements */
            tensor(&zpar, iscom, jscom, lscom, ismap, npair, alpha); /* Spin-dependent part (Clebsch-Gordan etc.) */

            if (sznz < 0) /* If SzNz type needs finalization from idpars */
              sznz = sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom); /* Reset N in ixcom/jxcom if changed */

            symnsq(nsqj, ins, iscom, jscom, &zpar); /* N(N+1) and N.S dependence */
            parskp = (fabs(zpar) < 1e-30); /* If K-independent part is effectively zero */
            if (parskp) continue; /* Skip this parameter */

            if (si2 < 0) { /* Commutator with N*N/2, for certain distortion terms */
                           /* This part might be ([H, N^2]/2N_i)*P_i = (N_f^2-N_i^2)/2N_i * P_i * <f|op|i>, but is likely simpler */
                           /* Effective operator is multiplied by (N(N+1)_bra - N(N+1)_ket)/2, but here it's simplified or handled by rmatrx/symnsq? */
                           /* For si2<0 (commutator), rmatrx returns 1. symnsq calculates (N(N+1)_bra - N(N+1)_ket)/2 effectively. */
                           /* The following line seems to imply it becomes N_bra * Op, which is unusual for a commutator form. Needs clarification. */
                           /* Pickett's J. Mol. Spec. 148, 371-377 (1991) Eq. A5 implies such operators are related to N_alpha for specific TYP. */
                           /* Here si2 = -1 comes from idpari if specific type of operator. Let's assume zpar already has the core operator part. */
                           /* The logic in symnsq seems to handle the N(N+1) difference for inq > 0. */
                           /* If si2 < 0, it means the operator form is like {A, N^2} or [A, N^2]. For this specific case (si2<0), it simplifies. */
                           /* This specific multiplication by 'ni' seems to be an ad-hoc scaling or part of a specific operator definition. */
                           /* From documentation: if si2 is negative, it implies a commutator with N*N for some parameter types. */
                           /* The multiplication by 'ni' may be specific to how 'zpar' is constructed for these commutator-type operators. */
              zpar *= ni;
            }
          } /* End of 'else' for new type of K-dependent operator or non-MNSQ */

          /* TODO change here: */
          if (ni == 0 && nsqj > 0 && ins == 0 && sznz == 0) continue; /* N(N+1) term for N=0 state is zero if no spins changing N effectively */

          ele = zpar * sqj[nsqj] * spar_now->zfac; /* Combine K-indep (zpar), N(N+1)^nsqj (sqj[nsqj]), and overall scaling (spar_now->zfac) */
          if (isgn < 0) /* Apply overall sign factor */
            ele = -ele;

          i = ipder[ipar]; ipbase = ipar; /* Get derivative index for this parameter */
          if (i < 0) { /* If parameter is constrained (i < 0 points to main parameter) */
            ipbase = -1 - i; /* Index of the parameter it's constrained to */
            i = ipder[ipbase]; /* Derivative index of the main parameter */
            ele *= par[ipar]; /* Multiply by the value of the constrained parameter (ratio) */
          }

          if (first) { /* First pass: Hamiltonian construction */
            pbar += ele * par[ipbase]; /* Add contribution to the composite parameter pbar */
            oldpar = TRUE; /* Current K-dependent part might be reusable */
          } else { /* Second pass: Derivatives */
            /* Add contribution to derivative: dE/dP_main = Sum_states <psi| (dH_op / dP_main) |psi> */
            /* dbar[k] contains <psi_k| (K-dep_op * K-indep_op_coeffs_excluding_P_main) |psi_k> */
            /* So, dedp for P_main gets ele * dbar (where ele has P_constr / P_main if constrained) */
            daxpy(nsize, ele, dbar, 1, &dedp[i * ndm], 1); /* dedp[i*ndm + k] += ele * dbar[k] for each state k */
          }
        } /* end parameter loop (spar_now) */
      } /* end JXX loop (ket sub-blocks) */
    } /* end IXX loop (bra sub-blocks) */

    if (first) { /* End of first pass (Hamiltonian construction) */
      if (ifdump) { /* If only dumping the Hamiltonian matrix (no diagonalization) */
        for (i = 0; i < nsize; ++i) {
          pt = &t[i + i * ndm]; /* Diagonal element */
          *pt += egy0;          /* Add E0 offset */
          egy[i] = (*pt);       /* Store diagonal element as 'energy' */
          pmix[i] = 0.;         /* No mixing */
        }
        return 0; /* Successfully dumped */
      }
      if (glob.idiag < 0 || nsize <= 1) { /* If no diagonalization requested or trivial 1x1 block */
        nd = nsize + 1;
        dcopy(nsize, t, nd, egy, 1); /* Copy diagonal elements of H to egy */
        for (i = 0; i < nsize; ++i) { /* Create identity eigenvector matrix */
          pt = &t[i * ndm];
          dcopy(nsize, &zero, 0, pt, 1); /* Zero out column */
          pmix[i] = pt[i] = 1.;          /* Eigenvector is unit vector, pmix=1 */
        }
        if (glob.idiag == 4 && nsize == 1) /* Special case for K-sort of 1x1 */
          pmix[0] = (double) kbgni; /* pmix stores K value */
      } else { /* Diagonalization is requested and nsize > 1 */
        if (glob.idiag == 0) { /* Standard energy sort, check for K-roll before diagonalization */
          roll = kroll(nsize, t, nsblk, ibkptr, ikmin); /* Modifies t if roll-over found */
        } else if (glob.idiag == 2 || glob.idiag == 5) { /* Sort by input H diagonal order */
          nd = nsize + 1;
          dcopy (nsize, t, nd, dbar, 1); /* Save original diagonal of H in dbar for later sorting */
        }
        /* diagonalize Hamiltonian matrix t, eigenvalues in egy, eigenvectors in t (overwriting H) */
#if HAM_DEBUG /* Debug block to check matrix before/after diagonalization */
        nd = nsize * nsize;
        pscr = (double *) mallocq((size_t)nd * sizeof(double)); /* Scratch space */
        if (pscr == NULL && nd > 0) { perror("mallocq failed for pscr"); exit(EXIT_FAILURE); }
        if (nd > 0) pscr[0] = t[0]; /* To use pscr if nd > 0 */
        dcopy(nd, t, 1, pscr, 1); /* Save original H */
#endif
        nqnsep = hdiag(nsize, nsize, t, egy, pmix, iqnsep); /* Diagonalize */
#if HAM_DEBUG /* Debug block to check T E T_dagger = H_diag and T T_dagger = I */
        if (pscr != NULL) { /* Check pscr was allocated */
          ele = 0.;
          for (i = 0; i < nsize; ++i) {
            for (ii = 0; ii <= i; ++ii) {
              pbar = (nd > 0) ? -pscr[i + ii * nsize] : 0.0; /* Original H_ii element */
              for (ixx = 0; ixx < nsize; ++ixx) { /* Sum E_k * T_ik * T_jk */
                jxx = ixx * nsize;
                pbar += egy[ixx] * t[i + jxx] * t[ii + jxx];
              }
              ele += fabs(pbar); /* Accumulate difference | H_orig_ij - (T E T_dagger)_ij | */
              pbar = ddot(nsize, &t[i] ,nsize, &t[ii], nsize); /* Check T_i . T_j */
              if (i == ii) pbar -= 1.0; /* Should be 0 for orthogonality */
              ele += fabs(pbar); /* Accumulate orthogonality error */
            }
          }
          free(pscr); pscr = NULL;
        }
#endif
        if (nqnsep < 0) { /* Error during diagonalization */
          printf(" diagonalization failure, err = %d in block %d, %d\n",
                 nqnsep, iblk, nsize);
          exit(EXIT_FAILURE);
        }
      	switch (glob.idiag) { /* Sorting options for eigenvalues/vectors */
      	default:  /* idiag = 0: energy sort of subblock (Wang sub-blocks are sorted internally) */
      	  ordblk(nsize, nsize, iqnsep, t, egy, ibkptr, pmix, idx); /* idx is scratch */
	        break;
        case 1:  /* projection sort of full block (sort by <K> or similar from pmix) */
          for (i = 0; i <= nsize; ++i) { /* Create identity permutation for full block sort */
            jdx[i] = (short) i;
          }
          ordblk(nsize, nsize, iqnsep, t, egy, jdx, pmix, idx); /* idx is scratch */
          break;
        case 2: /* Energy sort within Wang sub-blocks, but ensure it follows original H diagonal order if energies are close */
          i = ordblk(nsize, nsize, iqnsep, t, egy, ibkptr, pmix, idx); /* Initial sort by energy within sub-blocks */
          if (ODD2(i)) { /* ODD2(i) might mean some degeneracy or near-degeneracy was handled by ordblk */
            ordham(nsize, iqnsep, dbar, ibkptr, jdx); /* Get permutation jdx based on original H diagonal (dbar) */
            fixham(nsize, nsize, t, egy, pmix, jdx); /* Apply this permutation */
          }
          break;
        case 3: /* Sort by tau = Ka-Kc like proxy, within vibrational and spin sub-blocks */
        case 4: /* Sort by K (projection of K^2), within vibrational and spin sub-blocks */
        case 5: /* Sort by original H diagonal, within vibrational and spin sub-blocks */
          vbkptr = &ibkptr[glob.maxblk]; /* Use area after ibkptr for vibrational sub-block pointers */
          nd = 0;
          ivbase = -1; /* Current vibrational state index being processed */
          for (i = 0; i < nsblk; ++i) { /* Determine boundaries of vibrational sub-blocks */
            ixx = ibkptr[i]; /* Start of current (Wang*spin) sub-block */
            /* Group sub-blocks by common vibrational state (and l-doublet component) */
            jvbase = (int) ((unsigned) ivs[i] >> 2) - (itau[i] & 1); /* Effective vib index for grouping l-doublets */
            if (ivbase != jvbase) { /* If new vibrational group starts */
              vbkptr[nd++] = (short) ixx; /* Store start index of this new group */
              ivbase = jvbase;
            }
          }
          vbkptr[nd] = (short) nsize; /* Mark end of last vibrational group */
          i = ordblk(nsize, nsize, iqnsep, t, egy, vbkptr, pmix, idx); /* Initial energy sort within these vib groups */
          if ((i & 2) == 0 /* TODO INCLUDE?: && glob.idiag != 4 */) /* If no issues from ordblk, and not K-sort, we might be done */
            break;                            /* For idiag=4 (bestk), always proceed */
          if (glob.idiag == 5){ /* Sort by original H diagonal within vib groups */
            ordham(nsize, iqnsep, dbar, vbkptr, idx); /* Get permutation based on H_diag (dbar) within vib groups */
            fixham(nsize, nsize, t, egy, pmix, idx); /* Apply permutation */
            break;
          }
          /* For idiag=3 (tau) or idiag=4 (K) */
          for (i = 0; i < nsblk; ++i) { /* Fill idx with sorting keys (tau values stored in itau earlier) */
            ixx = itau[i];
            iz = ibkptr[i + 1];
            for (ii = ibkptr[i]; ii < iz; ++ii) {
              idx[ii] = (short) ixx; /* Assign sorting key for each basis state */
              ixx += 8; /* Increment key for next basis state in sub-block, arbitrary large step */
            }
          }
          if (glob.idiag == 4) { /* Sort by K values */
            /* bestk calculates K expectation values into wk, gets permutation jdx, applies to t, egy, pmix */
            bestk(nsize, nsize, iqnsep, vbkptr, idx, jdx, t, egy, pmix, wk);
            dcopy(nsize, wk, 1, pmix, 1); /* Store K expectation values in pmix */
            fixham(nsize, nsize, t, egy, pmix, jdx); /* Apply final permutation for K sort */
            break;
          }
          /* idiag == 3 (tau sort) */
          for (i = 0; i < nsize; ++i) { /* Copy tau keys from idx to wk (double) */
            wk[i] = (double) idx[i];
          }
          ordham(nsize, iqnsep, wk, vbkptr, idx); /* Get permutation based on tau values (in wk) */
          fixham(nsize, nsize, t, egy, pmix, idx); /* Apply permutation */
          break;
        } /* end switch (glob.idiag) */
      } /* end else (diagonalization requested) */

      if (isneg) { /* Correct energies for oblate representation by negating */
      	for (i = 0; i < nsize; ++i)
	        egy[i] = -egy[i];
      }
      for (i = 0; i < nsize; ++i) /* Add back the E0 offset to all energies */
	      egy[i] += egy0;
      first = FALSE; /* Mark first pass (Hamiltonian calculation and diagonalization) as done */
    } else { /* Not first pass (first == FALSE), means this is the derivative calculation pass */
      if (roll) { /* If K-roll occurred, derivatives might be unreliable. This seems to recalculate energies from derivatives? */
                  /* This part seems to be a way to get energies if K-roll modified H, by summing P * (dE/dP)_approx */
                  /* However, dE/dP are usually calculated FROM final E and T. This might be a fallback or alternative. */
                  /* More likely, this sums P * (d<H>/dP) if eigenvectors are identity (no diagonalization case). */
                  /* If roll=TRUE, it means KROLL modified the H matrix. The derivatives d(H_modified)/dP are being used. */
                  /* This will effectively give the expectation value of H_modified if par are the parameters. */
      	dcopy(nsize, &zero, 0, egy, 1); /* Zero out energy array */
	      for (ipar = npar - 1; ipar >= 0; --ipar) { /* Loop over all parameters */
	        i = ipder[ipar]; /* Get derivative index */
	        if (i >= 0) { /* If not a constrained parameter */
	          ele = par[ipar]; /* Parameter value */
	          daxpy(nsize, ele, &dedp[i * ndm], 1, egy, 1); /* egy[k] += par[ipar] * dedp(k, ipar) */
                                                            /* This accumulates Sum_p (P * dE/dP) into egy. If dE/dP are correct first-order derivatives, this would be E. */
          }
        }
      }
      break; /* Exit the 'do..while(!first)' loop after the second pass (derivatives) */
    }
  } while (!first);  /* This loop runs twice: first=TRUE, then first=FALSE */

  i = 1;
  if (idpar[0] == (bcd_t)0) /* Check if parameter list is empty (idpar[0] stores length) */
    i = -1;     /* This is just to ensure idpar is "used" to avoid compiler warnings. Return value is effectively status. */
  return i; /* Return status */
}  /* hamx */

/**
 * @brief Clears a 2D array (or a block of it) by setting elements to zero.
 * @param n1 Number of rows (or elements if ix=0 or vec is 1D conceptually).
 * @param n2 Number of columns (if ix=1, effectively n1*n2 total elements).
 * @param vec Pointer to the array/vector to be cleared.
 * @param ix Stride/increment: 0 for contiguous, 1 for column-major like access if n1 is leading dim.
 *           (dcopy behavior: if ix=0, all n1 elements are zeroed. if ix=1, elements vec[0], vec[1]... zeroed.
 *            If this is meant for a matrix (N1 rows, N2 cols) with ix=1 for column-major,
 *            it would clear the first N1*N2 elements contiguously.
 *            However, dcopy(count, src, inc_src, dest, inc_dest) means this is dcopy(N, &zero, 0, pvec, ix).
 *            If ix=1, it clears N elements contiguously. If ix=0, it's not standard for dcopy.
 *            Given dcopy(&zero,0,...), src value is always zero, inc_src is irrelevant.
 *            This function clears n1*n2 elements of vec, advancing by ix for each of n1*n2 elements.
 *            This seems to be intended to clear an n1*n2 matrix when ix=1 means contiguous.
 *            If ix > 1, it would clear strided.
 *            The usage in hamx is dclr(nsize, nsize, t, 1) and dclr(nsize, n, dedp, 1),
 *            where 1 means contiguous clearing of nsize*nsize or nsize*n elements.
 *            This dclr itself calls dcopy(count, &zero, 0, pvec, ix).
 *            If ix=1 in dclr, then dcopy uses incdest=1. If ix=0 in dclr, dcopy uses incdest=0 (setting only pvec[0]).
 *            The hamx calls use ix=1, so it's contiguous clearing.
 *            The 'ix' parameter in dclr seems to be directly passed as incDest to dcopy.
 *            So dclr(N, M, vec, 1) clears N*M elements of vec contiguously.
 *            dclr(N, M, vec, stride) clears N*M elements with specified stride.
 *            Actually, dcopy takes (count, src_ptr, src_inc, dest_ptr, dest_inc).
 *            Here, dcopy((int)nsq, &zero, 0, pvec, ix) will set (int)nsq elements of pvec (starting at pvec[0])
 *            to zero, where each element is pvec[k*ix]. This seems to be clearing with a stride 'ix'.
 *            In hamx, ix=1 is used, meaning contiguous clearing.
 * )
 * @return int Always 0.
 */
int dclr(n1, n2, vec, ix)
const int n1, n2, ix;
double *vec;
{    /*  clear a N1*N2 block */ /* Original comment */
  static long nbig_chunk_size = 0x7ff0; /* Max elements for dcopy at once (approx 32k, for 16-bit int limit?) */
  long n_elements_to_clear; /* Renamed nsq */
  double *current_vec_ptr;  /* Renamed pvec */

  current_vec_ptr = vec;
  n_elements_to_clear = n1 * (long) n2; /* Total number of elements to clear */
  while (n_elements_to_clear > nbig_chunk_size) /* If total exceeds chunk size for dcopy */
  {
    dcopy((int)nbig_chunk_size, &zero, 0, current_vec_ptr, ix);
    current_vec_ptr += nbig_chunk_size; /* <--- Ensure this simple advance is used */
    n_elements_to_clear -= nbig_chunk_size;
  }
  dcopy((int)n_elements_to_clear, &zero, 0, current_vec_ptr, ix); /* Clear remaining elements */
  return 0;
} /* dclr */

/**
 * @brief Applies Euler series transformation to operator matrix elements.
 * @param neuler Type of Euler operator (determines 'a' or 'b' denominator and power). Odd values set aden/bden, even values apply them.
 * @param newblk Pointer to a flag, set to TRUE if vibrational states change (to re-init aden/bden).
 * @param nsqj Pointer to power of N(N+1). Modified by this function.
 * @param ikq Pointer to power of K^2. Modified by this function.
 * @param ksi K value for bra side of the sub-block start.
 * @param ksj K value for ket side of the sub-block start.
 * @param ni N quantum number for bra state.
 * @param nj N quantum number for ket state.
 * @param ncos Number of elements in wk, idx, jx.
 * @param wk Array of K-dependent matrix element parts (modified in place).
 * @param ix Array of row indices (relative to sub-block start).
 * @param jx Array of column indices (relative to sub-block start).
 * @param par Parameter value (aden or bden if neuler is odd).
 * @return int Number of non-zero elements after transformation (ncos), or 0 if transformation results in zero.
 */
int specop(neuler, newblk, nsqj, ikq, ksi, ksj, ni, nj, ncos, wk, ix, jx, par)
const int neuler, ksi, ksj, ni, nj, ncos;
int *nsqj, *ikq;
BOOL *newblk;
double *wk;
const double par;
const short *ix, *jx;
{
#define NSPOP 5 /* Number of Euler series denominators (a_n, b_n up to n=4, since ist = (neuler-2)/2) */
  static double aden[NSPOP], bden[NSPOP]; /* Denominators 'a' and 'b' for Euler series */
  double akisq, akjsq, anisq, anjsq, da, dd, di, dj, xi, xj; /* K^2_bra, K^2_ket, (N(N+1)-K^2)_bra, (N(N+1)-K^2)_ket, a-d, d, denominator_bra, denominator_ket, element_bra_part, element_ket_part */
  double sqi, sqj, di0, dj0; /* N(N+1)_bra, N(N+1)_ket, (1+d*N(N+1))_bra, (1+d*N(N+1))_ket */
  int i_elem, nkq, nkq0, nnq, nnq0, ki, kj, ist; /* i_elem: loop counter; nkq,nkq0: K^2 power current,original; nnq,nnq0: N(N+1) power current,original; ki,kj: K_bra,K_ket for current element; ist: Euler series index (0 to NSPOP-1) */
  unsigned int ipwr; /* Loop counter for applying powers */

  if (*newblk) { /* If new vibrational block, reset denominators */
    for (i_elem = 0; i_elem < NSPOP; ++i_elem) {
      aden[i_elem] = 0.;
      bden[i_elem] = 0.;
    }
    *newblk = FALSE;
  }
  ist = (neuler - 2) >> 1; /* Determine which Euler series (0 for type 2,3; 1 for 4,5 etc.) */
  i_elem = ist;
  if (i_elem >= NSPOP) i_elem = 0; /* Safety, though neuler should ensure this */
  nkq0 = (*ikq);  /* Original K^2 power from parameter ID */
  nnq0 = (*nsqj); /* Original N(N+1) power from parameter ID */

  if (ODD(neuler)) { /* Odd neuler types are used to set the denominator values */
    if (nnq0 == 0 && nkq0 == 1) { /* If parameter is for K^2 term of Euler series (e.g., a_k * K^2 / (1 + a_k K^2 + b_k N(N+1))) -> sets 'a' */
      aden[i_elem] = par;
    } else if (nnq0 == 1 && nkq0 == 0) { /* If parameter is for N(N+1) term of Euler series -> sets 'b' */
      bden[i_elem] = par;
    }
    return 0; /* No modification to wk, just setting denominators */
  }
  /* Even neuler types apply the transformation using stored aden/bden */
  dd = bden[i_elem]; da = aden[i_elem] - dd; /* d = b_k; a-d = a_k - b_k */
  sqi = (double)(ni * (ni + 1)); di0 = 1. + dd * sqi; /* Denom_bra_const_part = 1 + b_k * N_bra(N_bra+1) */
  if (ni == nj) {
    sqj = sqi; dj0 = di0;
  } else {
    sqj = (double)(nj * (nj + 1)); dj0 = 1. + dd * sqj; /* Denom_ket_const_part */
  }

  for (i_elem = 0; i_elem < ncos; ++i_elem) { /* Loop over non-zero K-dependent elements */
    nkq = nkq0; /* Reset current K^2 power */
    nnq = nnq0; /* Reset current N(N+1) power */
    xi = wk[i_elem]; /* Current K-dependent part of matrix element */
    ki = ksi + ((int) ix[i_elem] << 1); /* K_bra for this specific element */
    akisq = (double) (ki * ki);         /* K_bra^2 */
    di = di0 + da * akisq;              /* Full denominator for bra side: 1 + b_k N(N+1) + (a_k-b_k)K^2 = 1 + a_k K^2 + b_k (N(N+1)-K^2) */
    if (fabs(di) < 1e-9) di = (di > 0.0) ? 1e-9 : -1e-9; /* Avoid division by zero */
    akisq /= di; /* K_bra^2 / Denom_bra (effective K^2 term) */

    kj = ksj + ((int) jx[i_elem] << 1); /* K_ket for this specific element */
    if (ki == kj && ni == nj) { /* Diagonal in K and N */
      if (nnq > 0) {            /* (N(N+1)-K^2) part */
        ipwr = (unsigned int) nnq;
        anisq = sqi / di - akisq; /* (N_bra(N_bra+1)/Denom_bra) - (K_bra^2/Denom_bra) */
        for (;;) {
          if (ODD(ipwr)) xi *= anisq;
          if ((ipwr >>= 1) == 0) break;
          anisq *= anisq;
        }
      }
      if (nkq > 0) {            /* K^2 part */
        ipwr = (unsigned int) nkq;
        for (;;) {
          if (ODD(ipwr)) xi *= akisq; /* akisq is already K^2/D here */
          if ((ipwr >>= 1) == 0) break;
          akisq *= akisq;
        }
      }
      if (ist != 0) /* For higher Euler terms (neuler >= 4), an extra 1/D factor */
        xi /= di;
      wk[i_elem] = xi;
    } else { /* Off-diagonal in K or N */
      akjsq = (double) (kj * kj); /* K_ket^2 */
      dj = dj0 + da * akjsq;      /* Full denominator for ket side */
      if (fabs(dj) < 1e-9) dj = (dj > 0.0) ? 1e-9 : -1e-9;
      akjsq /= dj; /* K_ket^2 / Denom_ket (effective K^2 term) */
      xj = xi; /* Ket part starts with same base K-dependent element */

      if (nnq > 0) {            /* (N(N+1)-K^2) part */
        ipwr = (unsigned int) nnq;
        anisq = sqi / di - akisq; /* Bra side (N(N+1)-K^2)/D */
        anjsq = sqj / dj - akjsq; /* Ket side (N(N+1)-K^2)/D */
        for (;;) {
          if (ODD(ipwr)) {
            xi *= anisq;
            xj *= anjsq;
          }
          if ((ipwr >>= 1) == 0) break;
          anisq *= anisq;
          anjsq *= anjsq;
        }
      }
      if (nkq > 0) {            /* K^2 part */
        ipwr = (unsigned int) nkq;
        for (;;) {
          if (ODD(ipwr)) {
            xi *= akisq; /* akisq is K_bra^2/D_bra */
            xj *= akjsq; /* akjsq is K_ket^2/D_ket */
          }
          if ((ipwr >>= 1) == 0) break;
          akisq *= akisq;
          akjsq *= akjsq;
        }
      }
      if (ist != 0) { /* Extra 1/D factor for higher Euler terms */
        xi /= di;
        xj /= dj;
      }
      wk[i_elem] = 0.5 * (xi + xj); /* Average of bra and ket transformed parts */
    }
  }
  *nsqj = 0; /* Powers are now incorporated into wk */
  *ikq = 0;  /* Powers are now incorporated into wk */
  return ncos;
}                               /* specop */

/**
 * @brief Applies Fourier coefficient scaling to operator matrix elements.
 * @param ifc Fourier coefficient type/index. ifc=0 initializes/stores rho. Negative ifc implies sine series.
 * @param iv Vibrational index for bra state.
 * @param jv Vibrational index for ket state.
 * @param kdel Change in K for the operator.
 * @param ksi K value for bra side of the sub-block start.
 * @param ksj K value for ket side of the sub-block start.
 * @param ncos Number of elements in wk, idx, jx.
 * @param wk Array of K-dependent matrix element parts (modified in place).
 * @param ix Array of row indices (relative to sub-block start).
 * @param jx Array of column indices (relative to sub-block start).
 * @return int Always 0.
 */
int specfc(ifc, iv, jv, kdel, ksi, ksj, ncos, wk, ix, jx)
const int ifc, iv, jv, kdel, ksi, ksj, ncos;
double *wk;
const short *ix, *jx;
{
  static double pfac[MAXVIB];  /* Stores rho_v * pi/3 for each vibrational state v */
  static short ipfac[MAXVIB]; /* Stores sign flag for rho_v (1 if original rho_v was negative) */
  double di_rho_term, dj_rho_term, angle_arg; /* Renamed di, dj, ang */
  int i_elem, fourier_order_n, ki, kj, kd_elem, is_g12_term_active, k_min_for_g12; /* Renamed i, n, kmin */

  if (ifc == 0) { /* Initialization call or storing rho */
    if (ncos == 0) { /* Global initialization call (e.g. from pasort) */
      for (i_elem = 0; i_elem < MAXVIB; ++i_elem) {
        pfac[i_elem] = 0.;
        ipfac[i_elem] = 0;
      }
    } else { /* Storing rho_v for vibrational state 'iv' (wk[0] contains rho_v) */
      /*  multiply by pi/3 */
      pfac[iv] = wk[0] * acos(0.5); /* acos(0.5) = pi/3 */
      if (pfac[iv] < 0.) { /* Store absolute value and sign separately */
        pfac[iv] = -pfac[iv];
        ipfac[iv] = 1; /* Flag that original rho_v was negative */
      }
    }
    return 0;
  }
  fourier_order_n = ifc;
  if (fourier_order_n < 0) /* Negative ifc indicates sine series */
    fourier_order_n = -1 - fourier_order_n; /* Map to positive index for series order, e.g. -1 -> 0, -2 -> 1 */

  is_g12_term_active = glob.g12 & fourier_order_n; /* Check if G12 symmetry term (-1)^K is active for this Fourier order */
                                               /* glob.g12 is likely 1 if G12 option active, 0 otherwise. */
                                               /* This makes is_g12_term_active = 1 if G12 active AND fourier_order_n is odd, else 0. */
                                               /* However, the doc says "if K is odd", not if fourier_order_n is odd. This might be a simplified check. */
                                               /* If G12=1, then for odd fourier_order_n, the (-1)^K term is applied. */

  angle_arg = (double) fourier_order_n; /* Base Fourier order n for cos(n*angle) or sin(n*angle) */
  /* Adjust base angle if rho_v and rho_v' have different signs (ipfac) -> corresponds to cos/sin((n-1/2)*angle) */
  if (ODD(ipfac[iv] + ipfac[jv])) angle_arg -= 0.5;

  di_rho_term = angle_arg * pfac[iv]; /* n * rho_v' * pi/3 */
  dj_rho_term = angle_arg * pfac[jv]; /* n * rho_v'' * pi/3 */

  for (i_elem = 0; i_elem < ncos; ++i_elem) {
    ki = ksi + (ix[i_elem] << 1); /* K_bra for current element */
    kj = ksj + (jx[i_elem] << 1); /* K_ket for current element */
    kd_elem = ki - kj; k_min_for_g12 = kj; /* k_min_for_g12 is K of lower state if K changes, or just K if K doesn't change. */
    if (kd_elem < 0) {
      kd_elem = -kd_elem; k_min_for_g12 = ki;
    }
    if (kdel != kd_elem) /* This condition seems to handle parity of K-functions or phase, reverses sign of kj for angle calc */
      kj = -kj;         /* if the operator's Delta K doesn't match element's Delta K (e.g. for P+ vs P- type components) */
                      /* This effectively changes the relative phase of the K-dependent angle for one side. */

    if (fourier_order_n != 0) { /* If not a constant term (n=0) */
      if ((is_g12_term_active & k_min_for_g12) != 0) /* Apply G12 term: if G12 active for this fourier_order_n AND k_min_for_g12 is odd, negate wk[i] */
                                                /* This interpretation of (is_g12_term_active & k_min_for_g12) seems to imply is_g12_term_active is a bitmask rather than just 0/1. */
                                                /* If glob.g12=1, is_g12_term_active is 1 if fourier_order_n is odd. So this is (-1)^K_min if fourier_order_n is odd. */
        wk[i_elem] = -wk[i_elem];

      angle_arg = di_rho_term * ki + dj_rho_term * kj; /* Argument for cos/sin: (n*rho'_K'_eff + n*rho''_K''_eff) * pi/3, where K_eff can be +/-K */
                                                   /* If documentation implies K_avg, this form (K_bra_term + K_ket_term) is equivalent if rho_bra = rho_ket */
                                                   /* Or it's (n * rho_v' * K_v' + n * rho_v'' * K_v'')*pi/3 form directly */
      if (ifc < 0) { /* Sine series */
        wk[i_elem] *= sin(angle_arg);
      } else { /* Cosine series */
        wk[i_elem] *= cos(angle_arg);
      }
    } else { /* For ifc=0 (or fourier_order_n = 0), it's a derivative term like K_avg * operator */
             /* From IDPAR docs: ifc=0, NS=0, TYP=91, KSQ=0, V1=V2 gives rho. */
             /* This ifc here is from IDPAR's FF field. if FF=0, it's usually no Fourier term. */
             /* An FF value of 0 (for cosine) or 10 (for sine) might lead to fourier_order_n=0 here. */
             /* Doc implies TYP=91 with KSQ=0, V1=V2 gives rho. Parameter ID for that would have FF=0. */
             /* The parameter is "rho". If this routine is called with ifc leading to fourier_order_n=0, it's likely applying a K-dependent factor. */
             /* The operator is multiplied by (K_bra + K_ket)/2. */
      wk[i_elem] *= 0.5 * (ki + kj);
    }
  }
  return 0;
} /* specfc */

/**
 * @brief Fixes/adjusts N quantum numbers in ixcom/jxcom and spin values in iscom/jscom
 *        when SzNz type operators effectively change N.
 * @param sznz Flag indicating type of SzNz operator and which N to modify.
 *             Values 1,2,3 modify bra (ixcom, iscom). Values 4,5 modify ket (jxcom, jscom).
 *             Negative values indicate reset after operation.
 * @param ni N quantum number for bra state (may be modified if sznz affects bra).
 * @param nj N quantum number for ket state (may be modified if sznz affects ket).
 * @param ixcom Quantum number block for bra state (XNVAL may be modified).
 * @param jxcom Quantum number block for ket state (XNVAL may be modified).
 * @param iscom Spin quantum numbers for bra state (iscom[nspin] i.e. 2N may be modified).
 * @param jscom Spin quantum numbers for ket state (jscom[nspin] i.e. 2N may be modified).
 * @return int 0 if N values were reset, or -sznz if N values were changed for an operation.
 *         Returns 0 (and does nothing) if the operation is invalid (e.g., N becomes <0).
 */
int sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom)
const int sznz; /* sznz_operator_type_flag */
const int ni, nj; /* N_bra, N_ket (original, before modification by this call) */
int *ixcom, *jxcom, *iscom, *jscom;
{
  if (sznz <= 0) { /* Negative or zero sznz means reset N to original values */
    if (sznz <= -4) { /* Reset ket side (nj) */
      jscom[nspin] = nj << 1; /* jscom[nspin] stores 2*N_ket */
      jxcom[XNVAL] = nj;
    } else { /* Reset bra side (ni) or both if sznz == 0 (no-op but resets state) */
      iscom[nspin] = ni << 1; /* iscom[nspin] stores 2*N_bra */
      ixcom[XNVAL] = ni;
    }
    return 0; /* Return 0 to indicate reset */
  } else { /* Positive sznz means apply operator that changes N */
    switch (sznz) {
    case 1: /* S_z N_z, no change to N, but check if N=0 valid */
      if (ni + nj == 0 && ni == 0) /* If both N_bra and N_ket are 0 initially, this operator is likely zero unless spins allow change */
        return 0; /* Invalid if N=0 and no mechanism to change it */
      break;
    case 2: /* Operator like N_-, changes N_bra to N_bra-1 */
      if (ni <= 0) /* N_bra cannot be 0 or less to decrease */ /* Original check ni <=1, N=0 means N_z is usually 0 */
        return 0;
      iscom[nspin] -= 2; /* iscom[nspin] is 2*N_bra, so decrement by 2*(1) */
      --ixcom[XNVAL];    /* Decrement N_bra in quantum block */
      break;
    case 3: /* Operator like N_+, changes N_bra to N_bra+1 */
      if (ni < 0) /* N_bra cannot be negative */ /* Original check ni <=0, but N=0 can be raised */
        return 0;
      iscom[nspin] += 2;
      ++ixcom[XNVAL];
      break;
    case 4: /* Operator like N_-, changes N_ket to N_ket-1 (from ket perspective) */
      if (nj <= 0) /* N_ket cannot be 0 or less to decrease */ /* Original check nj <=1 */
        return 0;
      jscom[nspin] -= 2;
      --jxcom[XNVAL];
      break;
    case 5: /* Operator like N_+, changes N_ket to N_ket+1 (from ket perspective) */
      if (nj < 0) /* N_ket cannot be negative */ /* Original check nj <= 0 */
        return 0;
      jscom[nspin] += 2;
      ++jxcom[XNVAL];
      break;
    default: /* Unknown sznz type */
      return 0;
    }
    return -sznz; /* Return negative of original sznz to indicate modification occurred */
  }
} /* sznzfix */

/**
 * @brief Applies SzNz operator scaling to K-dependent matrix elements.
 * @param ni N quantum number for bra state (potentially modified by sznzfix).
 * @param nj N quantum number for ket state (potentially modified by sznzfix).
 * @param ksi K value for bra side of the sub-block start.
 * @param ksj K value for ket side of the sub-block start.
 * @param iscom Spin quantum numbers for bra state. iscom[nspin] is 2*N_bra. iscom[0] is 2*J_bra. iscom[nspin+1] is 2*S_bra.
 * @param jscom Spin quantum numbers for ket state. jscom[nspin] is 2*N_ket. jscom[0] is 2*J_ket. jscom[nspin+1] is 2*S_ket.
 * @param ncos Number of elements in wk, idx, jx.
 * @param wk Array of K-dependent matrix element parts (modified in place).
 * @param ix Array of row indices (relative to sub-block start).
 * @param jx Array of column indices (relative to sub-block start).
 * @return int Always 0.
 */
int sznzop(ni, nj, ksi, ksj, iscom, jscom, ncos, wk, ix, jx)
const int ni, nj, ksi, ksj, ncos;
const int *iscom, *jscom;
double *wk;
const short *ix, *jx;
{
  int iflg_n_change_type, i_elem, k_val, kk_val, n_bra_times_2, n_ket_times_2, n_bra_orig_times_2, n_ket_orig_times_2, j_qn_times_2, spin_qn_times_2; /* Renamed iflg, i, k, kk, nni, nnj, nni0, nnj0, jj, iis */
  double val_bra, val_ket, current_wk_element, sum_contrib; /* Renamed vali, valj, twk, sum */

  if (ncos <= 0) /* No elements to operate on */
    return 0;

  n_bra_times_2 = ni << 1; /* Current 2*N_bra (possibly modified by sznzfix) */
  n_bra_orig_times_2 = iscom[nspin]; /* Original 2*N_bra from spin table */
  val_bra = (double) (n_bra_orig_times_2 + 1); /* Factor related to (2N_orig+1) */

  n_ket_times_2 = nj << 1; /* Current 2*N_ket */
  n_ket_orig_times_2 = jscom[nspin]; /* Original 2*N_ket */
  val_ket = (double) (n_ket_orig_times_2 + 1); /* Factor related to (2N_orig_ket+1) */

  if (n_bra_times_2 != n_bra_orig_times_2) { /* If N_bra was changed by sznzfix (e.g. N -> N-1) */
    iflg_n_change_type = 1; /* Flag: N_bra changed */
    val_bra = sqrt(val_bra * (n_bra_times_2 + 1)); /* Scaling factor like sqrt((2N_orig+1)(2N_final+1)) for <N_final|N_op|N_orig> */
  } else if (n_ket_times_2 != n_ket_orig_times_2) { /* If N_ket was changed */
    iflg_n_change_type = -1; /* Flag: N_ket changed */
    val_ket = sqrt(val_ket * (n_ket_times_2 + 1));
  } else { /* N_bra and N_ket were not changed by sznzfix (e.g. pure SzNz, no N ladder operator) */
    iflg_n_change_type = 0;
  }

  /* Calculate ket-side Wigner coefficient part (if N_ket is involved or N_bra is not changed) */
  if (iflg_n_change_type <= 0) {
    j_qn_times_2 = jscom[0];         /* 2*J_ket */
    spin_qn_times_2 = jscom[nspin + 1]; /* 2*S_ket (first spin) */
    /* Phase factor related to K_ket and angular momenta values */
    k_val = ksj + ksj + n_ket_orig_times_2 + n_ket_times_2 + j_qn_times_2 + spin_qn_times_2;
    if ((k_val & 3) != 0) /* Check parity of sum of (2K_ket + 2N_orig_ket + 2N_final_ket + 2J_ket + 2S_ket), effectively a phase */
      val_ket = -val_ket;
    /* Multiply by 6j symbol: (N_orig_ket S_ket J_ket / N_final_ket S_ket 1) where 1 is tensor order of N_z */
    val_ket *= c6jj(n_ket_orig_times_2, 2, n_ket_times_2, spin_qn_times_2, j_qn_times_2, spin_qn_times_2);
  }
  /* Calculate bra-side Wigner coefficient part (if N_bra is involved or N_ket is not changed) */
  if (iflg_n_change_type >= 0) {
    j_qn_times_2 = iscom[0];         /* 2*J_bra */
    spin_qn_times_2 = iscom[nspin + 1]; /* 2*S_bra */
    k_val = ksi + ksi + n_bra_orig_times_2 + n_bra_times_2 + j_qn_times_2 + spin_qn_times_2;
    if ((k_val & 3) != 0)
      val_bra = -val_bra;
    val_bra *= c6jj(n_bra_times_2, 2, n_bra_orig_times_2, spin_qn_times_2, j_qn_times_2, spin_qn_times_2); /* Note order of N's for bra */
  }

  for (i_elem = 0; i_elem < ncos; ++i_elem) { /* Loop over K-dependent elements */
    current_wk_element = wk[i_elem];
    sum_contrib = 0.;
    if (iflg_n_change_type <= 0) { /* Contribution from ket side (or both if iflg=0) */
      k_val = ksj + ((int) jx[i_elem] << 1); /* K_ket for current element */
      kk_val = k_val << 1;                   /* 2*K_ket */
      /* Multiply by K_ket and 3j symbol: <N_final_ket, K_ket | N_z | N_orig_ket, K_ket> ~ K_ket * (N_final_ket 1 N_orig_ket / -K_ket 0 K_ket) */
      sum_contrib = val_ket * current_wk_element * k_val * c3jj(n_ket_orig_times_2, 2, n_ket_times_2, -kk_val, 0, kk_val);
    }
    if (iflg_n_change_type >= 0) { /* Contribution from bra side (or add to ket if iflg=0) */
      k_val = ksi + ((int) ix[i_elem] << 1); /* K_bra for current element */
      kk_val = k_val << 1;                   /* 2*K_bra */
      sum_contrib += val_bra * current_wk_element * k_val * c3jj(n_bra_times_2, 2, n_bra_orig_times_2, -kk_val, 0, kk_val);
    }
    wk[i_elem] = sum_contrib; /* Update element with SzNz scaling */
  }
  return 0;
} /* sznzop */

/**
 * @brief Calculates the overall symmetry of an interaction block.
 * @param ixcom Quantum number block for the bra state.
 * @param jxcom Quantum number block for the ket state.
 * @return unsigned int The symmetry of the interaction (0=A, 1=Bx, 2=By, 3=Bz in D2),
 *                      result of XORing bra and ket state symmetries.
 */
unsigned int blksym(ixcom, jxcom)
const int *ixcom, *jxcom;
{                               /* get block symmetry from state symmetry */ /* Original comment */
  return (unsigned int) (3 & (ixcom[XSYM] ^ jxcom[XSYM])); /* XOR symmetries and mask to 2 bits */
} /* blksym */

/**
 * @brief Orders eigenvalues within sub-blocks to roughly match the diagonal of the Hamiltonian.
 *        Used when glob.idiag = 2 or 5 for initial ordering before fixham.
 * @param nn Dimension of the Hamiltonian.
 * @param mask Eigenvector mixing coefficient or projection quantum number (not directly used for sorting here, but passed).
 * @param egy Array of eigenvalues (sorted by energy within sub-blocks initially by hdiag),
 *            also contains original diagonal elements of H if idiag=2,5. This function sorts 'egy' based on its values.
 * @param isblk Array of pointers to the start of each sub-block (Wang or vibrational).
 * @param iswap Output array storing the permutation required to achieve the desired order.
 * @return int Always 0.
 */
int ordham(nn, mask, egy, isblk, iswap)
const int nn;
double *egy; /* Contains original diagonal H elements if idiag=2,5; energies if idiag=3,4 (where wk is passed as egy) */
const short *isblk; /* Pointers to start of sub-blocks */
short *iswap, *mask; /* mask is iqnsep from hamx, used here to check mcmp */
{
/* subroutine to order eigenvalues within sub-block like diagonal of */ /* Original comment */
/*      Hamiltonian */ /* Original comment */
  double temp_egy_val; /* Renamed tmp */
  int current_sub_block_idx, i_row, i_check, i_min_idx, sub_block_end_idx, inext_isblk_val; /* Renamed iblk, i, ii, iq, is, inext */
  short mcmp_val; /* Renamed mcmp */

  current_sub_block_idx = 1; /* sub-block counter, seems 1-based from isblk[1] */
  i_row = 0; /* Current row/state index we are trying to place */
  inext_isblk_val = isblk[1]; /* Index of the start of the *next* sub-block */
  for (sub_block_end_idx = 1; sub_block_end_idx < nn; ++sub_block_end_idx) { /* Loop through states, up to nn-1 */
    if (sub_block_end_idx == inext_isblk_val) { /* If we reached the end of the current sub-block */
      inext_isblk_val = isblk[++current_sub_block_idx]; /* Get end of next sub-block */
      iswap[i_row] = (short) i_row; /* The last element of a sub-block is already in place relative to itself */
    } else {                    /* Within a sub-block, find element with minimum 'egy' value */
      temp_egy_val = egy[i_row]; /* Current minimum 'egy' value */
      i_min_idx = i_row;         /* Index of that minimum */
      mcmp_val = mask[i_row];    /* Get mask (iqnsep) value for current state */
      for (i_check = sub_block_end_idx; i_check < inext_isblk_val; ++i_check) { /* Search remaining elements in this sub-block */
        if (egy[i_check] < temp_egy_val && mask[i_check] == mcmp_val) { /* If smaller 'egy' found and same mask type */
          temp_egy_val = egy[i_check];
          i_min_idx = i_check;
        }
      }
      iswap[i_row] = (short) i_min_idx; /* Record where the element for row 'i_row' was found */
      if (i_min_idx > i_row) { /* If the minimum was not at 'i_row', swap it into 'i_row' */
        egy[i_min_idx] = egy[i_row];
        egy[i_row] = temp_egy_val;
      }
    }
    i_row = sub_block_end_idx; /* Move to next row to process */
  }
  /* The last element (nn-1) needs its iswap value set */
  if (nn > 0) iswap[nn-1] = (short)(nn-1); /* Assuming it's correctly placed by previous logic or only one left */

  return 0;
} /* ordham */

/**
 * @brief Reorders eigenvalues, eigenvectors, and mixing coefficients based on a permutation vector.
 * @param ndm Leading dimension of matrix t (eigenvectors).
 * @param nn Dimension of the Hamiltonian / number of states.
 * @param t Eigenvector matrix (modified in place).
 * @param egy Eigenvalue array (modified in place).
 * @param p Mixing coefficient array (pmix) (modified in place).
 * @param iswap Permutation array from ordham or bestk. iswap[i] is the original index of the state that should go to position i.
 * @return int Always 0.
 */
int fixham(ndm, nn, t, egy, p, iswap)
const int ndm, nn;
double *t, *egy, *p;
const short *iswap;
{
/* subroutine to order eigenvalues within sub-block like diagonal of */ /* Original comment (slightly misleading, it applies a pre-calculated order) */
/*      Hamiltonian using permutation found with ordham */ /* Original comment */
  int i_current, original_idx_for_current_pos; /* Renamed i, iq */
  int loop_end_val; /* Renamed is */

  loop_end_val = nn - 2; /* Loop from nn-2 down to 0 */
  if (loop_end_val >= 0) {
    for (i_current = loop_end_val; i_current >= 0; --i_current) {
      original_idx_for_current_pos = iswap[i_current]; /* Get original index of the element that should be at i_current */
      if (original_idx_for_current_pos > i_current) { /* If it's not already in place, swap it */
        /* etswap swaps eigenvectors, eigenvalues, and pmix values between original_idx_for_current_pos and i_current */
        etswap(ndm, nn, original_idx_for_current_pos, i_current, t, egy, p);
      }
    }
  }
  return 0;
}   /* fixham */

/**
 * @brief Performs K-based sorting of eigenvalues and eigenvectors, calculating K expectation values.
 * @param ndm Leading dimension of matrix t.
 * @param nsize Dimension of the Hamiltonian.
 * @param iqnsep Eigenvector mixing coefficients or projection labels from hdiag.
 * @param ibkptr Array of pointers to the start of each vibrational sub-block.
 * @param itau Array of tau-like quantum numbers (Ka-Kc proxy) or K values for each basis state.
 * @param iswap Output array for the permutation vector.
 * @param t Eigenvector matrix (modified in place).
 * @param egy Eigenvalue array (modified in place).
 * @param pmix Mixing coefficient array (modified in place, will store <K_op^2> / <sum |c_i|^2>).
 * @param wk Work array.
 * @return int Always 0.
 */
int bestk(ndm, nsize, iqnsep, ibkptr, itau, iswap, t, egy, pmix, wk)
    const int ndm,
    nsize;
short *iqnsep, *ibkptr, *itau, *iswap;
double *t, *egy, *pmix, *wk;
{
#define MAXNK 4
  double ele;
  long ndml, tindx, tindx0;
  int i, iz, kd, k, ii, jj, ibgn, iend, iblk, is, iq, nk;
  short mcmp;
  ndml = ndm;
  dcopy(nsize, &zero, 0, wk, 1);
  for (i = 0; (ibgn = ibkptr[i]) < nsize; ++i)
  {
    iend = ibkptr[i + 1] - 1;
    tindx0 = ibgn * (ndml + 1);
    for (ii = ibgn; ii <= iend; ++ii)
    {
      /* calculate expectation for PzPz */
      kd = (int)(itau[ii] >> 2);
      tindx = tindx0;
      for (jj = ibgn; jj <= iend; ++jj)
      {
        ele = kd * t[tindx];
        wk[jj] += ele * ele;
        tindx += ndml;
      }
      ++tindx0;
    }
  }
  for (i = 0; i < nsize; ++i)
  {
    wk[i] = sqrt(wk[i] / pmix[i]);
  }
  iend = ibkptr[1];
  iblk = 1;
  i = 0;
  for (is = 1; is < nsize; ++is)
  {
    if (is == iend)
    {
      iend = ibkptr[++iblk];
    }
    else
    { /* find min wk and itau */
      mcmp = iqnsep[i];
      iq = i;
      ele = wk[i];
      iz = i;
      jj = itau[i];
      for (ii = is; ii < iend; ++ii)
      {
        if (mcmp != iqnsep[ii])
          continue;
        if (wk[ii] < ele)
        {
          ele = wk[ii];
          iq = ii;
        }
        if (itau[ii] < (short)jj)
        {
          iz = ii;
          jj = itau[ii];
        }
      }
      if (iq > i)
      {
        wk[iq] = wk[i];
        wk[iq] = ele;
        etswap(ndm, nsize, iq, i, t, egy, pmix);
      }
      iswap[i] = (short)iz;
      if (iz > i)
      {
        itau[iz] = itau[i];
        itau[i] = (short)jj;
      }
    }
    i = is;
  }
  iend = ibkptr[1];
  iblk = 1;
  i = 0;
  nk = 0;
  for (is = 1; is < nsize; ++is)
  {
    if (is == iend)
    {
      iend = ibkptr[++iblk];
      nk = 0;
    }
    else
    {
      if (nk <= 0)
      {
        kd = itau[i] >> 2;
        nk = 1;
        for (ii = is; ii < iend; ++ii)
        {
          k = itau[ii] >> 2;
          if (k > kd)
            break;
          /* find basis with K == kd, nk = the number of basis with same K */
          ++nk;
        } /* loop ii */
      }
      if (nk > 1)
      {
        mcmp = iqnsep[i];
        ele = egy[i];
        iq = i;
        for (iz = 1; iz < nk; ++iz)
        {
          ii = iz + i;
          if (egy[ii] < ele && mcmp == iqnsep[ii])
          {
            ele = egy[ii];
            iq = ii;
          }
        }
        if (iq > i)
        {
          ele = wk[iq];
          wk[iq] = wk[i];
          wk[iq] = ele;
          etswap(ndm, nsize, iq, i, t, egy, pmix);
        }
      }
      --nk;
    }
    i = is;
  }
  return 0;
} /* bestk */

/**
 * @brief Checks for Hamiltonian roll-over (diagonal elements not monotonically increasing/decreasing).
 *        If roll-over is detected, it modifies the Hamiltonian to enforce monotonicity
 *        and zeros out corresponding off-diagonal elements.
 * @param hamiltonian_dim Leading dimension of matrix t (Hamiltonian). Renamed nsizd.
 * @param t Hamiltonian matrix (modified in place if roll-over).
 * @param num_sub_blocks Number of sub-blocks. Renamed nsblk.
 * @param sub_block_pointers Array of pointers to the start of each sub-block. Renamed sbkptr.
 * @param kmin_values Array of minimum K values for each sub-block. Renamed kmin.
 * @return BOOL TRUE if roll-over was detected and corrected, FALSE otherwise.
 */
BOOL kroll(nsizd, t, nsblk, sbkptr, kmin)
    const int nsizd,
    nsblk;
double *t;
const short *sbkptr, *kmin;
{
  /* subroutine to make sure that diagonal elements of the Hamiltonian */
  /*    are monotonically increasing (decreasing for oblate) */
  /*    after K=KTHRSH */
  double vall, tlast, tmp, *ptmp;
  long ndmt;
  int ibgn, iend, i, k, n, ixx;
  BOOL roll;

  ndmt = nsizd + 1;
  roll = FALSE;
  for (ixx = 0; ixx < nsblk; ++ixx)
  {
    i = 5 - kmin[ixx];
    ibgn = sbkptr[ixx];
    if (i >= 2)
      ibgn += i >> 1;
    iend = sbkptr[ixx + 1] - 1;
    if (ibgn >= iend)
      continue;
    ptmp = &t[ibgn * ndmt];
    tmp = tlast = *ptmp;
    ++ibgn;
    for (i = ibgn; i <= iend; ++i)
    {
      ptmp += ndmt;
      tlast = tmp;
      tmp = (*ptmp);
      if (tmp < tlast)
        break;
    }
    if (i > iend)
      continue;
    k = kmin[ixx] + ((i - sbkptr[ixx]) << 1);
    printf(" roll-over at K = %d\n", k);
    roll = TRUE;
    vall = 1e15;
    for (; i <= iend; ++i)
    {
      n = i - 1;
      dcopy(n, &zero, 0, &t[i], nsizd);
      n = nsizd - i;
      dcopy(n, &zero, 0, &t[i * ndmt], 1);
      tlast += vall;
      *ptmp = tlast;
      ptmp += ndmt;
    }
  }
  return roll;
} /* kroll */

/**
 * @brief Gets quantum numbers for a given sub-block basis state.
 * @param packed_vib_spin_sym Packed (Vibrational_index, Symmetry, Spin_pattern_index). Renamed im.
 * @param f_qn_times_2 2*F quantum number for the block. Renamed iff.
 * @param sub_block_dim Dimension of this specific sub-block. Renamed nsiz.
 * @param k_begin Starting K value for this sub-block. Renamed kbgn.
 * @param ixcom_out Output array (NDXCOM size) to store rotational quanta. Renamed ixcom.
 * @param iscom_out Output array (MAXNS size) to store coupled angular momenta (2N, 2J, 2F1...). Renamed iscom.
 * @param vib_idx_out Output vibrational state index. Renamed iv.
 * @return int Spin pattern index (ispx) for this state.
 */
int getqs(im, iff, nsiz, kbgn, ixcom, iscom, iv)
    const int im,
    iff, nsiz, kbgn;
int *ixcom, *iscom, *iv;
{
  /*     subroutine to get quanta corresponding to sub-block */
  /*     on entry: */
  /*         IM= position in sub-block mold */
  /*         IFF = F quantum number *2 */
  /*         NSIZ= size of sub-block */
  /*         KBGN= first K value */
  /*     on return: */
  /*         IXCOM = list of rotational containing */
  /*                 0: DIMENSION OF SUB-BLOCK */
  /*                 1: SYMMETRY CODE FOR BLOCK (0=A,1=Bx,2=By,3=Bz) */
  /*                 2: 'N' QUANTUM NUMBER */
  /*                 3: BEGINNING K */
  /*                 4: L VALUE */
  /*                 5: VIBRATIONAL QUANTA */
  /*                 5: ITOT SYMMETRY */
  /*                 5: EXTRA ITOT QUANTA */
  /*         ISCOM = list of angular momenta */
  /*         IV = vibrational quantum number   if L=0,1 */
  /*         IV = vibrational quantum number-1 if L=-1 */
  /*     find list of angular momenta */
  SVIB *pvinfo;
  int *isscom;
  short *iis, *jjs;
  int i, iiv, nset, nspinv, ispx, k, lupper;
  unsigned int mvs;

  mvs = (unsigned int)im;
  ispx = (int)(mvs >> glob.msshft);
  iiv = (int)(mvs >> 2) & glob.msmask;
  pvinfo = &vinfo[iiv];
  iis = pvinfo->spt;
  nset = iis[0];
  jjs = &iis[nset * ispx];
  i = (int)(mvs & 3);
  iscom[0] = iff + (*jjs);
  ixcom[XNVAL] = iscom[0] >> 1; /* find N */
  ixcom[XDIM] = nsiz;
  ixcom[XSYM] = i; /* find symmetry */
  if (nsiz < 0)
    setgsym((int)pvinfo->gsym);
  k = kbgn;
  if (k < 0)
  {
    k = pvinfo->knmin[i];
    i = pvinfo->knmin[3 - i];
    if (k > i)
      k = i;
    if (k < 0)
      k = 0;
  }
  ixcom[XKBGN] = k;
  ixcom[XVIB] = iiv;
  ixcom[XISYM] = 0;
  lupper = pvinfo->lvupper;
  if (ODD(lupper))
    --iiv;
  *iv = iiv;
  i = pvinfo->lvqn;
  ixcom[XLVAL] = i;
  if (i != 0 || k != 0)
    i = 1;
  ixcom[XIQN] = i;
  if (nspin == 0)
    return 1;
  nspinv = nset - 1;
  if (nspinv > nspin)
    nspinv = nspin;
  isscom = &iscom[nspin];
  isscom[0] = iscom[0];
  for (i = 1; i < nspinv; ++i)
  {
    iscom[i - 1] = iff + jjs[i];
    isscom[i] = iis[i];
  }
  iscom[i - 1] = iff;
  isscom[i] = iis[i];
  while (i < nspin)
  {
    ++i;
    iscom[i - 1] = iff;
    isscom[i] = 0;
  }
  if (glob.nitot > 0)
  {
    if (glob.nitot >= 3)
    {
      if (itsym < nspinv)
      {
        i = jjs[itsym + 1];
        iscom[itsym] = i;
        k = MOD(i >> 2, glob.nitot);
        i -= k << 2;
        if (is_esym[k] != 0)
          ++i;
        ixcom[XISYM] = k;
        ixcom[XIQN] |= i;
      }
      else
      {
        iscom[itsym] = 0;
      }
    }
    iscom[itptr] = jjs[itptr + 1];
  }
  return ispx;
} /* getqs */

/**
 * @brief Parses a parameter identifier (from SPAR struct) into its constituent fields.
 * @param pspar Pointer to the SPAR structure containing the packed parameter ID and pre-parsed fields.
 * @param ksq_out Output: power of K^2.
 * @param itp_out Output: Euler series type/index (pspar->euler).
 * @param tensor_order_N_out Output: tensor order in N. Renamed ln.
 * @param dir_cos_order_out Output: order of direction cosine part. Renamed ld.
 * @param k_delta_out Output: change in K. Renamed kdel.
 * @param n_dot_s_power_out Output: power of N.S. Renamed ins.
 * @param spin1_idx_out Output: index of first spin. Renamed si1.
 * @param spin2_idx_out Output: index of second spin. Renamed si2.
 * @param sznz_type_out Output: SzNz operator type. Renamed sznz.
 * @param fourier_coeff_idx_out Output: Fourier coefficient index. Renamed ifc.
 * @param alpha_sym_out Output: I_tot alpha symmetry component. Renamed alpha.
 * @param l_delta_out Output: change in l. Renamed ldel.
 * @param k_avg_out Output: K_avg value if applicable. Renamed kavg.
 * @return int Flags from the pspar structure.
 */
int idpars(pspar, ksq_out, itp_out, tensor_order_N_out, dir_cos_order_out, k_delta_out, n_dot_s_power_out, spin1_idx_out, spin2_idx_out, sznz_type_out,
           fourier_coeff_idx_out, alpha_sym_out, l_delta_out, k_avg_out) /* Renamed parameters for clarity */
SPAR *pspar;
int *ksq_out, *itp_out, *tensor_order_N_out, *dir_cos_order_out, *k_delta_out, *n_dot_s_power_out, *spin1_idx_out, *spin2_idx_out, *sznz_type_out; /* Renamed parameters */
int *fourier_coeff_idx_out, *alpha_sym_out, *l_delta_out, *k_avg_out; /* Renamed parameters */
{
/*     subroutine to parse parameter identifiers to subfields */ /* Original comment */
/*     on entry: */ /* Original comment */
/*         IDPAR= 10*( parameter ID/1000) + ITP */ /* Original comment - this refers to an older BCD packing, now fields are directly in pspar */
/*     on return: */ /* Original comment */
/*         KSQ= power of K*K */ /* Original comment */
/*         ITP= parameter subtype */ /* Original comment - now pspar->euler is primary type, itp_out is that */
/*         LN= tensor order */ /* Original comment */
/*         LD= order of direction cosine ( L <= LD ) */ /* Original comment */
/*         KDEL= change in K */ /* Original comment */
/*         INS=power of N.S */ /* Original comment */
/*         SI1,SI2=spin identifiers */ /* Original comment */
  int flags_return; /* Renamed iret */
  *ksq_out = (int) pspar->ksq;
  *itp_out = (int) pspar->euler; /* This is the primary Euler type/flag now */
  *tensor_order_N_out = (int) pspar->mln;
  *dir_cos_order_out = (int) pspar->mld;
  *k_delta_out = (int) pspar->mkdel;
  *n_dot_s_power_out = (int) pspar->mins;
  *spin1_idx_out = (int) pspar->msi1;
  *spin2_idx_out = (int) pspar->msi2;
  *sznz_type_out = (int) pspar->msznz;
  *fourier_coeff_idx_out = (int) pspar->fc;
  *alpha_sym_out = (int)pspar->alpha;
  *l_delta_out = (int)pspar->mldel;
  *k_avg_out = (int)pspar->kavg;
  flags_return = (int)pspar->flags;
  return flags_return;
} /* idpars */

/**
 * @brief Determines tensor orders for spin couplings and checks triangle inequalities.
 * @param llf_total_tensor_order 2 * Tensor order for the total F quantum number.
 * @param dir_cos_order Order of the direction cosine operator part.
 * @param n_tensor_order Order of the N-dependent tensor part.
 * @param k_delta_op Change in K for the operator.
 * @param spin1_idx Index of the first spin involved.
 * @param spin2_idx Index of the second spin involved.
 * @param lscom_tensor_orders_out Output array: lscom_out[0] to lscom_out[nspin-1] are 2*tensor orders for intermediate couplings (L_J, L_F1, ...).
 *                  lscom_out[nspin] to lscom_out[nspin+nspin] are 2*tensor orders for individual spin operators (L_S, L_I1, ...).
 * @param iscom_bra_qns Array of 2*angular momentum values for the bra state.
 * @param jscom_ket_qns Array of 2*angular momentum values for the ket state.
 * @return int Maximum spin index involved in the operator (up to nspin), or -1 if any triangle rule is violated.
 */
int getll(llf_total_tensor_order, dir_cos_order, n_tensor_order, k_delta_op, spin1_idx, spin2_idx, lscom_tensor_orders_out, iscom_bra_qns, jscom_ket_qns)
const int llf_total_tensor_order,
dir_cos_order, n_tensor_order, k_delta_op, spin1_idx, spin2_idx;     /* Renamed llf, ld, ln, kd, si1, si2 */
int lscom_tensor_orders_out[];                                       /* Renamed lscom */
const int iscom_bra_qns[], jscom_ket_qns[];                          /* Renamed iscom, jscom */
{
  int *l_individual_spin_tensors_ptr;                                                                                                                                                                             /* Renamed lsscom */
  int i_spin_level, current_intermediate_L_times_2, ll_dir_cos_times_2, ll_n_tensor_times_2, max_L_N_operator_times_2, delta_2N_val, sum_2N_val, max_involved_spin_idx, error_return_flag, num_spins_minus_2_val; /* Renamed i, llj, lld, lln, llmax, jdif, jsum, maxspin, ierr, nm2 */

  error_return_flag = -1;
  max_L_N_operator_times_2 = 0;
  ll_dir_cos_times_2 = dir_cos_order << 1;
  ll_n_tensor_times_2 = n_tensor_order << 1;
  delta_2N_val = iscom_bra_qns[nspin] - jscom_ket_qns[nspin];
  if (delta_2N_val < 0)
    delta_2N_val = -delta_2N_val;

  if (ll_dir_cos_times_2 >= ll_n_tensor_times_2)
  {
    if (ll_n_tensor_times_2 < delta_2N_val)
      return error_return_flag;
    max_L_N_operator_times_2 = k_delta_op << 1;
    if (max_L_N_operator_times_2 < ll_dir_cos_times_2)
      max_L_N_operator_times_2 = ll_dir_cos_times_2;
  }
  else
  {
    if (ll_dir_cos_times_2 < delta_2N_val)
      return error_return_flag;
    max_L_N_operator_times_2 = k_delta_op << 1;
    if (max_L_N_operator_times_2 < ll_n_tensor_times_2)
      max_L_N_operator_times_2 = ll_n_tensor_times_2;
  }
  sum_2N_val = iscom_bra_qns[nspin] + jscom_ket_qns[nspin];
  if (max_L_N_operator_times_2 > sum_2N_val)
    return error_return_flag;

  if (nspin == 0)
  {
    lscom_tensor_orders_out[0] = ll_n_tensor_times_2;
    return 0;
  }

  max_involved_spin_idx = spin1_idx;
  if (spin1_idx > itsym && itsym < nspin)
    max_involved_spin_idx = nspin;

  l_individual_spin_tensors_ptr = &lscom_tensor_orders_out[nspin];
  l_individual_spin_tensors_ptr[0] = ll_n_tensor_times_2;
  for (i_spin_level = 1; i_spin_level <= nspin; ++i_spin_level)
  {
    l_individual_spin_tensors_ptr[i_spin_level] = 0;
  }
  if (spin1_idx > 0)
  {
    l_individual_spin_tensors_ptr[spin1_idx] = 2;
    if (spin2_idx > 0)
      l_individual_spin_tensors_ptr[spin2_idx] += 2;
  }

  current_intermediate_L_times_2 = l_individual_spin_tensors_ptr[0];
  num_spins_minus_2_val = nspin - 2;
  for (i_spin_level = 0; i_spin_level <= num_spins_minus_2_val; ++i_spin_level)
  {
    if (i_spin_level == itsym)
    {
      if (spin1_idx <= itsym)
      {
        current_intermediate_L_times_2 = 0;
      }
      else if (spin2_idx <= itsym || spin1_idx == spin2_idx)
      {
        current_intermediate_L_times_2 = l_individual_spin_tensors_ptr[spin1_idx];
      }
      i_spin_level = num_spins_minus_2_val;
    }
    else if (i_spin_level >= max_involved_spin_idx)
    {
      current_intermediate_L_times_2 = llf_total_tensor_order;
    }
    else
    {
      current_intermediate_L_times_2 -= l_individual_spin_tensors_ptr[i_spin_level + 1];
      if (current_intermediate_L_times_2 < 0)
        current_intermediate_L_times_2 = -current_intermediate_L_times_2;
    }
    lscom_tensor_orders_out[i_spin_level] = current_intermediate_L_times_2;
    delta_2N_val = iscom_bra_qns[i_spin_level] - jscom_ket_qns[i_spin_level];
    if (delta_2N_val < 0)
      delta_2N_val = -delta_2N_val;
    if (current_intermediate_L_times_2 < delta_2N_val)
      return error_return_flag;
    sum_2N_val = iscom_bra_qns[i_spin_level] + jscom_ket_qns[i_spin_level];
    if (current_intermediate_L_times_2 > sum_2N_val)
      return error_return_flag;
  }
  lscom_tensor_orders_out[nspin - 1] = llf_total_tensor_order;
  if (llf_total_tensor_order != 0)
    max_involved_spin_idx = nspin;

  return max_involved_spin_idx;
} /* getll */

/**
 * @brief Creates a mask for dircos to select allowed K and L transitions based on symmetry and operator type.
 * @param xbra_qns Quantum number block for the bra state. Renamed xbra.
 * @param xket_qns Quantum number block for the ket state. Renamed xket.
 * @param k_delta_op Total Delta K for the operator. Renamed kd.
 * @param l_delta_op Change in l quantum number (Delta l) for the operator. Renamed ldel.
 * @param operator_flags Operator property flags (MODD, MSYM2, MLZ, MIDEN from SPAR/DIRCOS). Renamed loff.
 * @param alpha_sym_component I_tot alpha symmetry component (0 if not I_tot or alpha=0 component). Renamed alpha.
 * @return int Mask value: bit 0 for (K'=K+DK, L'=L+DL), bit 1 for (K'=DK-K, L'=DL-L), bit 2 for (K'=K-DK, L'=L-DL).
 *         Returns 0 if no transitions are allowed.
 */
int getmask(const int *xbra, const int *xket, const int kd, const int ldel,
            const int loff, const int alpha)
{
  /* CREATE MASK FOR DIRCOS */
  /* KD IS TOTAL K QUANTIUM CHANGE */
  int ldif, lbra, lket, mask, kbra, kket, nitot, kk, ksbra, ksket;
  mask = 7;
  nitot = glob.nitot;
  if (nitot >= 3)
  { /* check Itot symmetry */
    ksbra = xbra[XISYM];
    ksket = xket[XISYM];
    kk = (xbra[XIQN] ^ xket[XIQN]) >> 1;
    if (TEST(loff, MIDEN))
    { /* Identity operator */
      if (kk != 0)
        return 0;
    }
    else if (EVEN(kk))
    { /* A1|A1, A1|E, E|E, A2|A2, ..*/
      if (TEST(loff, MSYM2))
        return 0;
    }
    else
    { /* kk is ODD */
      if (TEST(loff, MSYM2))
        return 0; /* <A2| I_alpha |E> */
      if (is_esym[ksbra] == 0 && is_esym[ksket] == 0)
        return 0; /*  <A1 |I_alpha| A2> == 0 */
    }
    mask = 0;
    kk = ksket - ksbra;
    if (alpha == 0)
    {
      if (kk == 0)
        mask = 5;
    }
    else if (alpha < nitot)
    {
      if (MOD(kk + alpha, nitot) == 0)
        mask = 1;
      if (MOD(kk - alpha, nitot) == 0)
        mask += 4;
    }
    else
    { /* quasi-diagonal */
      if (MOD(ksbra + ksket - alpha, nitot) == 0)
        mask = 2;
    }
    if (mask == 0)
      return mask;
  }
  /* mask bit 0 (val = 1): kbra = kket + kd, lbra = lket + ldel */
  /* mask bit 2 (val = 4): kbra = kket - kd, lbra = lket - ldel */
  /* mask bit 1 (val = 2): kbra = kd - kket, lbra = ldel - lket */
  lbra = xbra[XLVAL];
  lket = xket[XLVAL];
  ldif = lbra - lket;
  if (ldel == 0)
  { /* operator diagonal in l */
    if (ldif != 0)
      mask &= 2; /* clear bits 0 and 2 */
  }
  else
  { /* operator off-diagonal in l */
    if (ldif == 0)
    {
      mask &= 2; /* clear bits 0 and 2 */
    }
    else if (ldif != ldel)
    {
      mask &= 6; /* clear bit 0 */
    }
    else if (-ldif != ldel)
    {
      mask &= 3; /* clear bit 2 */
    }
  }
  if ((mask & 2) == 0)
    return mask; /* bit 1 not set */
  mask &= 5;     /* clear bit 1 temporarily */
  /* quasi-diagonal, KD = KBRA + KKET */
  if (lbra + lket != ldel)
    return mask;
  kbra = xbra[XKBGN];
  kket = xket[XKBGN];
  kk = kd - kbra - kket;
  if (kk < 0)
    return mask;
  if (EVEN(xket[XIQN]))
    kk -= 2;
  if (EVEN(xbra[XIQN]))
    kk -= 2;
  if (kk < 0)
    return mask;
  return (mask + 2);
} /* getmask */

/**
 * @brief Calculates N-dependent correction factors for reduced matrix elements.
 *        Used for operators where tensor order changes between N and overall.
 * @param dir_cos_order Tensor order of the direction cosine part. Renamed ld.
 * @param vib_tensor_order Effective tensor order of the vibrational part or target N-space operator. Renamed lv.
 * @param ixcom_bra Quantum number block for the bra state.
 * @param jxcom_ket Quantum number block for the ket state.
 * @return double The N-dependent correction factor. Returns 0 if selection rules violated.
 */
double rmatrx(dir_cos_order, vib_tensor_order, ixcom_bra, jxcom_ket) /* Renamed ld, lv, ixcom, jxcom */
const int dir_cos_order, vib_tensor_order;             /* Renamed ld, lv */
const int *ixcom_bra, *jxcom_ket; /* Renamed ixcom, jxcom */
{
  double delta_N_sq, n_sum_sq, N_val_as_double, n_sum_plus_1_as_double, factor;            /* Renamed dlsq, dnsq, dn, del, tmp */
  int n_sum_val, n_diff_val, N_bra_val, N_ket_val, min_tensor_order, current_tensor_order; /* Renamed nsum, ndif, ni, nj, lmin, lt */

  factor = 1.;
  if (dir_cos_order == vib_tensor_order) /* If tensor order doesn't change */
    return factor;

  N_bra_val = ixcom_bra[XNVAL];
  N_ket_val = jxcom_ket[XNVAL];
  n_sum_val = N_bra_val + N_ket_val + 1;
  n_sum_plus_1_as_double = (double)n_sum_val;                 /* N_bra + N_ket + 1 */
  n_sum_sq = n_sum_plus_1_as_double * n_sum_plus_1_as_double; /* (N_bra + N_ket + 1)^2 */

  if (dir_cos_order > vib_tensor_order)
  { /* Contracting tensor order (e.g. N-dependent part has lower order than dir.cos part) */
    min_tensor_order = vib_tensor_order;
    current_tensor_order = dir_cos_order;
    if (min_tensor_order == 0)
    { /* Special case for L_vib = 0 */
      min_tensor_order = 1;
      factor = N_bra_val * (double)(N_bra_val + 1) / n_sum_plus_1_as_double; /* N(N+1) / (N_bra+N_ket+1) factor */
    }
  }
  else
  { /* Expanding tensor order */
    min_tensor_order = dir_cos_order;
    current_tensor_order = vib_tensor_order;
    if (min_tensor_order == 0)
    { /* Special case for L_dir_cos = 0 */
      min_tensor_order = 1;
      factor = N_bra_val * (double)(N_bra_val + 1) * n_sum_plus_1_as_double; /* N(N+1) * (N_bra+N_ket+1) factor */
    }
  }
  if (current_tensor_order >= n_sum_val) /* If L_current > N_bra + N_ket (violates triangle rule for N_bra, N_ket, L_current) */
    return 0.;

  n_diff_val = N_bra_val - N_ket_val;
  /* n_sum_plus_1_as_double = (double) n_sum_val; -- already set */
  if (n_diff_val != 0)
    N_val_as_double = (double)(n_diff_val * n_diff_val); /* (N_bra - N_ket)^2 */
  else
    N_val_as_double = (double)(n_sum_val * n_sum_val); /* This was del = n_sum_plus_1_as_double if ndif=0, seems like error, should be 0 if ndif=0? */
                                                       /* Let's assume if ndif==0, the term involving del/dlsq later should vanish or be handled. */
                                                       /* For ndif=0, del is effectively (N_bra+N_ket+1)^2. */
                                                       /* This 'del' is actually only used if ndif != 0 in the loop. */

  while (current_tensor_order > min_tensor_order)
  {
    delta_N_sq = (double)(current_tensor_order * current_tensor_order); /* L_current^2 */
    factor *= 0.25 * (n_sum_sq - delta_N_sq);                           /* ( (N_bra+N_ket+1)^2 - L_current^2 ) / 4 */
    if (n_diff_val != 0)                                                /* Correction term if N_bra != N_ket */
      /* This term is 1 - (N_bra-N_ket)^2 / L_current^2. factor *= (1 - (N_bra-N_ket)^2 / L_current^2) */
      factor -= factor * N_val_as_double / delta_N_sq; /* Original: tmp -= tmp * del / dlsq; */
    --current_tensor_order;
  }
  return sqrt(factor); /* Returns sqrt of accumulated product */
} /* rmatrx */

/**
 * @brief Applies scaling for powers of N(N+1) and N.S operators.
 * @param n_sq_power Power of N(N+1). Renamed inq.
 * @param n_s_power Power of N.S. Renamed ins.
 * @param iscom_bra_qns Array of 2*angular momentum values for bra state.
 * @param jscom_ket_qns Array of 2*angular momentum values for ket state.
 * @param matrix_element_factor Pointer to the matrix element factor to be modified. Renamed zval.
 * @return int Always 0.
 */
int symnsq(n_sq_power, n_s_power, iscom_bra_qns, jscom_ket_qns, matrix_element_factor) /* Renamed parameters */
const int n_sq_power, n_s_power, *iscom_bra_qns, *jscom_ket_qns; /* Renamed inq, ins, iscom, jscom */
double *matrix_element_factor;                 /* Renamed zval */
{
  double n_n_plus_1_power_factor, dot_N_S_bra, dot_N_S_ket, factor_bra_side, factor_ket_side;          /* Renamed x, dotns, dotnsp, zr, zl */
  int delta_2N, first_spin_qn_bra_times_2, N_bra_times_2, N_ket_times_2, J_bra_times_2, J_ket_times_2; /* Renamed modf, ispn, nni, nnj, j */
  unsigned int current_power;                                                                          /* Renamed ipwr */

  N_bra_times_2 = iscom_bra_qns[nspin];
  N_ket_times_2 = jscom_ket_qns[nspin];
  delta_2N = N_bra_times_2 - N_ket_times_2; /* 2*(N_bra - N_ket) */
  if (n_sq_power <= 0)                      /* If no N(N+1) power, delta_2N is effectively 0 for N(N+1) part */
    delta_2N = 0;
  if (delta_2N == 0 && n_s_power == 0) /* No N(N+1) or N.S terms */
    return 0;

  factor_ket_side = factor_bra_side = 0.5 * (*matrix_element_factor); /* Split factor for bra and ket sides */

  if (delta_2N != 0)
  { /* N_bra != N_ket AND n_sq_power > 0. Operator is like (N(N+1))^k_bra - (N(N+1))^k_ket */
    /* This section scales the ket-side contribution if N_bra != N_ket */
    if (N_ket_times_2 == 0 && N_bra_times_2 != 0)
    {                       /* If N_ket = 0 (and N_bra != 0) */
      factor_ket_side = 0.; /* Ket side N(N+1) is 0 */
    }
    else if (N_bra_times_2 != 0 && N_ket_times_2 != 0)
    { /* Both non-zero, apply ratio */
      /* Scale by (N_ket(N_ket+1) / N_bra(N_bra+1)) ^ n_sq_power for ket side */
      n_n_plus_1_power_factor = N_ket_times_2 * (double)(N_ket_times_2 + 2) / (N_bra_times_2 * (double)(N_bra_times_2 + 2)); /* Ratio of N(N+1) terms */
      current_power = (unsigned int)n_sq_power;
      for (;;)
      {
        if (ODD(current_power))
          factor_ket_side *= n_n_plus_1_power_factor;
        if ((current_power >>= 1) == 0)
          break;
        n_n_plus_1_power_factor *= n_n_plus_1_power_factor;
      }
    }
    /* If N_bra == 0 and N_ket != 0, factor_bra_side should be 0 but is not explicitly set. It will be multiplied by N_bra(N_bra+1) later, effectively zeroing it. */
  }

  if (n_s_power > 0)
  {                                                                                                                                  /* If N.S term is present */
    J_bra_times_2 = iscom_bra_qns[0];                                                                                                /* 2*J_bra */
    first_spin_qn_bra_times_2 = iscom_bra_qns[nspin + 1];                                                                            /* 2*S_bra */
    n_n_plus_1_power_factor = (double)(first_spin_qn_bra_times_2 * (first_spin_qn_bra_times_2 + 2));                                 /* 4*S(S+1) for bra */
    dot_N_S_bra = 0.125 * ((J_bra_times_2 - N_bra_times_2) * (double)(J_bra_times_2 + N_bra_times_2 + 2) - n_n_plus_1_power_factor); /* (J(J+1) - N(N+1) - S(S+1))/2 for bra, times 1/4 from 2J etc. */
                                                                                                                                     /* Actual formula: (J(J+1) - N(N+1) -S(S+1))/2. Here terms are 2X, so ( (2J)(2J+2)/4 - (2N)(2N+2)/4 - (2S)(2S+2)/4 )/2 */
                                                                                                                                     /* = ( J(J+1) - N(N+1) - S(S+1) ) / 2. The 0.125 factor comes from (2X)(2X+2) = 4X(X+1). So (4J(J+1) - 4N(N+1) - 4S(S+1))/8. */

    J_ket_times_2 = jscom_ket_qns[0];                                                                /* 2*J_ket */
    first_spin_qn_bra_times_2 = jscom_ket_qns[nspin + 1];                                            /* 2*S_ket (reused variable name is for S_ket here) */
    n_n_plus_1_power_factor = (double)(first_spin_qn_bra_times_2 * (first_spin_qn_bra_times_2 + 2)); /* 4*S(S+1) for ket */
    dot_N_S_ket = 0.125 * ((J_ket_times_2 - N_ket_times_2) * (double)(J_ket_times_2 + N_ket_times_2 + 2) - n_n_plus_1_power_factor);

    current_power = (unsigned int)n_s_power;
    for (;;)
    {
      if (ODD(current_power))
      {
        factor_bra_side *= dot_N_S_bra;
        factor_ket_side *= dot_N_S_ket;
      }
      if ((current_power >>= 1) == 0)
        break;
      dot_N_S_bra *= dot_N_S_bra;
      dot_N_S_ket *= dot_N_S_ket;
    }
  }
  *matrix_element_factor = factor_bra_side + factor_ket_side; /* Sum contributions from bra and ket sides */
  return 0;
} /* symnsq */

/**
 * @brief Applies K^2 power scaling to K-dependent matrix elements.
 * @param k_sq_power Power of K^2. Renamed ikq.
 * @param K_bra_start K value for bra side of the sub-block start. Renamed ksi.
 * @param K_ket_start K value for ket side of the sub-block start. Renamed ksj.
 * @param num_elements Number of elements in work_array_elements, row_indices, col_indices. Renamed n.
 * @param work_array_elements Array of K-dependent matrix element parts (modified in place). Renamed wk.
 * @param row_indices Array of row indices (relative to sub-block start). Renamed ix.
 * @param col_indices Array of column indices (relative to sub-block start). Renamed jx.
 * @return int Number of non-zero elements (num_elements if any K values are non-zero, 0 otherwise if all K=0).
 */
int symksq(k_sq_power, K_bra_start, K_ket_start, num_elements, work_array_elements, row_indices, col_indices) /* Renamed parameters */
const int k_sq_power, K_bra_start, K_ket_start, num_elements; /* Renamed ikq, ksi, ksj, n */
double *work_array_elements;                /* Renamed wk */
short *row_indices, *col_indices;           /* Renamed ix, jx */
{
  double k_sq_bra_factor, k_sq_ket_factor, element_bra_part, element_ket_part; /* Renamed akisq, akjsq, xi, xj */
  int i_elem, K_bra_element, K_ket_element, non_zero_found_flag;               /* Renamed i, ki, kj, zflg */
  unsigned int current_k_sq_power;                                             /* Renamed nkq */

  non_zero_found_flag = 0;
  for (i_elem = 0; i_elem < num_elements; ++i_elem)
  {
    element_bra_part = work_array_elements[i_elem];
    K_bra_element = K_bra_start + ((int)row_indices[i_elem] << 1); /* K_bra for current matrix element */
    K_ket_element = K_ket_start + ((int)col_indices[i_elem] << 1); /* K_ket for current matrix element */
    current_k_sq_power = (unsigned int)k_sq_power;

    if (K_bra_element != K_ket_element)
    {                                      /* Off-diagonal in K */
      non_zero_found_flag = num_elements;  /* Mark that some K values are involved */
      element_ket_part = element_bra_part; /* Ket side starts with same base element */
      k_sq_bra_factor = (double)(K_bra_element * K_bra_element);
      k_sq_ket_factor = (double)(K_ket_element * K_ket_element);
      for (;;)
      {
        if (ODD(current_k_sq_power))
        {
          element_bra_part *= k_sq_bra_factor;
          element_ket_part *= k_sq_ket_factor;
        }
        if ((current_k_sq_power >>= 1) == 0)
          break;
        k_sq_bra_factor *= k_sq_bra_factor;
        k_sq_ket_factor *= k_sq_ket_factor;
      }
      work_array_elements[i_elem] = 0.5 * (element_bra_part + element_ket_part); /* Average for symmetrized K^n operator */
    }
    else if (K_bra_element != 0)
    { /* Diagonal in K, and K is not zero */
      non_zero_found_flag = num_elements;
      k_sq_bra_factor = (double)(K_bra_element * K_bra_element);
      for (;;)
      {
        if (ODD(current_k_sq_power))
          element_bra_part *= k_sq_bra_factor;
        if ((current_k_sq_power >>= 1) == 0)
          break;
        k_sq_bra_factor *= k_sq_bra_factor;
      }
      work_array_elements[i_elem] = element_bra_part;
    }
    else
    { /* K_bra = K_ket = 0, so K^n factor is 0 unless k_sq_power=0 */
      work_array_elements[i_elem] = (k_sq_power == 0) ? element_bra_part : 0.0;
      if (k_sq_power == 0 && fabs(element_bra_part) > 1e-30)
        non_zero_found_flag = num_elements;
    }
  }
  return non_zero_found_flag; /* Returns num_elements if any K dependent term was non-zero, else 0 */
} /* symksq */

/**
 * @brief Calculates the contribution of a single operator (defined by wk, idx, jdx)
 *        to the derivatives of eigenvalues with respect to a parameter.
 *        Essentially computes Sum_k T_ki * Op_ij * T_jk for each eigenstate (column of T).
 * @param matrix_size Dimension of the Hamiltonian.
 * @param derivative_accumulator_dp Output array to store d<H_op>/dP for each eigenstate. Renamed dp.
 * @param eigenvector_matrix Eigenvector matrix T. Renamed t.
 * @param num_op_elements Number of non-zero elements in the operator Op. Renamed n.
 * @param op_element_values Values of the non-zero elements of Op. Renamed wk.
 * @param op_row_indices Row indices for Op. Renamed idx.
 * @param op_col_indices Column indices for Op. Renamed jdx.
 * @param is_unit_op Flag: 1 if Op is a unit matrix type (diagonal with constant value, idx[0],jdx[0] define sub-block).
 *                   0 if Op is general sparse.
 * @return int Always 0.
 */
int dpmake(nsize, dp, t, n, wk, idx, jdx, isunit)
    const int nsize,
    n, isunit;
double *dp;
const double *t, *wk;
const short *idx, *jdx;
{ /*  find derivative contribution from sub-diagonal */
  const double *pt;
  double ele, tz;
  int i, k, iz, jz;

  if (isunit != 0)
  {
    iz = idx[0];
    jz = jdx[0];
    if (iz != jz)
    { /* off-diagonal unit matrix */
      for (k = 0; k < nsize; ++k)
      {
        dp[k] = 2. * ddot(n, &t[iz], 1, &t[jz], 1);
        iz += nsize;
        jz += nsize;
      }
    }
    else
    { /* diagonal unit matrix */
      for (k = 0; k < nsize; ++k)
      {
        pt = &t[iz];
        dp[k] = ddot(n, pt, 1, pt, 1);
        iz += nsize;
      }
    }
  }
  else
  {
    dcopy(nsize, &zero, 0, dp, 1);
    for (i = 0; i < n; ++i)
    {
      iz = idx[i];
      jz = jdx[i];
      ele = wk[i];
      if (iz != jz)
      {
        ele *= 2;
        for (k = 0; k < nsize; ++k)
        {
          dp[k] += t[iz] * t[jz] * ele;
          iz += nsize;
          jz += nsize;
        }
      }
      else
      {
        for (k = 0; k < nsize; ++k)
        {
          tz = t[iz];
          dp[k] += tz * tz * ele;
          iz += nsize;
        }
      }
    }
  }
  return 0;
} /* dpmake */

/**
 * @brief Calculates spin coupling factors (Wigner 6j, 9j symbols) for an operator.
 * @param matrix_element_factor_ptr Pointer to the matrix element factor to be modified by spin coeffs. Renamed zval.
 * @param iscom_bra_qns Array of 2*angular momentum values for bra state. Renamed iscom.
 * @param jscom_ket_qns Array of 2*angular momentum values for ket state. Renamed jscom.
 * @param lscom_tensor_orders Array of 2*tensor orders for intermediate and spin operators. Renamed lscom.
 * @param spin_index_map Maps operator spin index to actual spin index in iscom/jscom. Renamed imap.
 * @param num_coupled_pairs Number of spin coupling stages. Renamed npair.
 * @param alpha_sym_component I_tot alpha symmetry component. Renamed alpha.
 * @return int Always 0.
 */
int tensor(zval, iscom, jscom, lscom, imap, npair, alpha)
double *zval;
const int iscom[], jscom[], lscom[], imap[];
int npair, alpha;
{
  /*     subroutine to find spherical tensor spin coupling coefficients */
  /*     on entry: */
  /*         Z= initial multiplier */
  /*         NS= number of spins */
  /*         ISCOM,JSCOM= list of angular momenta *2 */
  /*         LLCOM= list of tensor orders *2 */
  /*         LSCOM = list of tensor orders for spins *2 */
  /*     on return: */
  /*         Z= modified multiplier */
  double zsq, z;
  int isgn, i, jjf, jji, nnf, llj, nni, lln, lls, ssf, ssi, ix, jj[9];

  if (npair <= 0)
    return 0;
  z = (*zval);
  isgn = 0;
  zsq = 1.;
  for (i = 0; i < npair; ++i)
  {
    ix = imap[i + i];
    lln = lscom[ix];
    nni = iscom[ix];
    nnf = jscom[ix];
    ix = imap[i + i + 1];
    lls = lscom[ix];
    ssi = iscom[ix];
    ssf = jscom[ix];
    if (alpha >= 0 && ix == itptr)
    {
      ix = itsym + nspin + 1;
      getzitot(&z, lls, iscom[ix], &lscom[ix],
               &iscom[itsym], &jscom[itsym], alpha, glob.nitot);
      i = nspin - 1;
    }
    llj = lscom[i];
    jji = iscom[i];
    jjf = jscom[i];
    if (lln != 0)
    {
      if (llj == 0)
      { /* T * U */
        isgn += jji + ssi + nnf;
        z *= c6jj(nni, lln, nnf, ssf, jji, ssi);
      }
      else if (lls == 0)
      { /* U = 1 */
        isgn += jjf + nni + ssi + llj;
        zsq *= jji + 1;
        zsq *= jjf + 1;
        z *= c6jj(nni, llj, nnf, jjf, ssi, jji);
      }
      else
      { /* T X U */
        isgn += 2;
        zsq *= 2.5;
        ix = (jjf + 1) * (llj + 1);
        zsq *= ix;
        zsq *= jji + 1;
        jj[0] = nnf;
        jj[1] = lln;
        jj[2] = nni;
        jj[3] = ssf;
        jj[4] = lls;
        jj[5] = ssi;
        jj[6] = jjf;
        jj[7] = llj;
        jj[8] = jji;
        z *= c9jj(jj);
      }
    }
    else if (lls != 0)
    { /* T = 1 */
      isgn += jji + nni + ssf + llj;
      zsq *= jji + 1;
      zsq *= jjf + 1;
      z *= c6jj(ssi, llj, ssf, jjf, nni, jji);
    }
  }
  if (ODD2(isgn))
    z = -z;
  *zval = z * sqrt(zsq);
  return 0;
} /* tensor */

/**
 * @brief Initializes dipole data from BCD encoded dipole identifiers.
 *        Parses each dipole ID, determines its type, symmetry, involved vibrational states,
 *        and spin operators. It also sets flags for imaginary dipoles based on operator
 *        order and selected phase conventions. Invalid or unsupported dipole types are flagged.
 * @param list_unit File pointer for writing warning messages.
 * @param ifdiag_out Output: Pointer to a BOOL that will be set to TRUE if the global diagonalization
 *                   option (glob.idiag) is non-negative, FALSE otherwise.
 * @param nsav_out Output: Pointer to an integer, always set to 1 in this function.
 * @param num_dipole_params Total number of dipole parameters to process.
 * @param bcd_dipole_ids Array of dipole identifiers in BCD format. The first element bcd_dipole_ids[0]
 *                       is expected to contain the number of BCD digits per ID.
 * @param is_imaginary_dipole_out Output array of integers, where each element is set to 1 if the
 *                                corresponding dipole is determined to be imaginary based on phase
 *                                conventions, 0 otherwise.
 * @return int Always returns 0 after processing.
 */
int setint(lu, ifdiag, nsav, ndip, idip, isimag)
FILE *lu;
BOOL *ifdiag;
const int ndip;
int *nsav, *isimag;
bcd_t *idip;
{ /*   subroutine to initialize intensity data */
  static char *bad_msg[] = {"\n", "(1) bad symmetry field\n",
                            "(2) vib out of range\n", "(3) unknown  dipole type\n",
                            "(4) transition between states of different weight\n",
                            "(5) alpha is not 0\n", "(6) bad spin field\n",
                            "(7) K is greater tnan Kmax\n"};
  SVIB *pvib1, *pvib2;
  short *iiv1, *iiv2;
  double zfac;
  unsigned int ijv;
  int i, icase, isym, ivdif, ii, iv1, iv2, lv1, lv2, ifac, lv, nsym;
  int si1, ibcd, ndecv, nbcd, ipty, iflg, kd, ld, ifc, imag, ldel;
  int nimag[16], ioff, kmin, nmin, nn, kavg, good_case, bad_dip;
  bcd_t itmp;
  *ifdiag = (glob.idiag >= 0);
  *nsav = 1;
  if (nddip > 0)
  {
    free(dipinfo);
    dipinfo = NULL;
  }
  dipinfo = (SDIP *)mallocq((size_t)ndip * sizeof(SDIP));
  dipinfo0.fac = 0.;
  memcpy(dipinfo, &dipinfo0, sizeof(SDIP));
  nddip = ndip;
  for (i = 0; i < 16; ++i)
  {
    nimag[i] = 0;
  }
  ifac = glob.vibfac + 1;
  ndecv = glob.vibdec;
  nbcd = (int)idip[0] & 0x7f;
  for (i = 0, ibcd = 0; i < ndip; ++i, ibcd += nbcd)
  {
    ii = (int)idip[ibcd + ndecv + 1];
    si1 = (ii >> 4) & 0x0f;
    ii &= 0x0f;
    if (ndecv > 2)
      ii = 100 * ii + bcd2i(idip[ibcd + 3]);
    if (ndecv > 1)
      ii = 100 * ii + bcd2i(idip[ibcd + 2]);
    itmp = idip[ibcd + 1];
    isym = (int)itmp & 0x0f;
    ii = 10 * ii + ((int)(itmp >> 4) & 0x0f);
    iv2 = ii / ifac;
    iv1 = ii - iv2 * ifac;
    ii = (int)idip[ibcd + ndecv + 2];
    icase = ii & 0x0f;
    ifc = (ii >> 4) & 0x0f;
    if (icase == 7 || icase == 8)
    {
      ifc = 10 * ifc + si1;
      si1 = 0;
    }
    else
    {
      icase += ifc * 10;
      ifc = 0;
    }
    kavg = 0;
    ivdif = iv1 - iv2;
    if (ivdif < 0)
    {
      ii = iv1;
      iv1 = iv2;
      iv2 = ii;
    }
    pvib1 = &vinfo[iv1];
    pvib2 = &vinfo[iv2];
    lv1 = pvib1->lvupper;
    lv2 = pvib2->lvupper;
    if (ODD(lv1))
      --iv1;
    if (ODD(lv2))
      --iv2;
    ipty = (lv1 ^ lv2) & 2;
    lv1 = pvib1->lvqn;
    lv2 = pvib2->lvqn;
    /* MAKE SURE VIBRATIONS ARE ORDERED */
    iflg = 1;
    bad_dip = 0;
    zfac = 1.;
    do
    { /* one pass */
      /*  CONVERT SYMMETRY */
      bad_dip = 1;
      if (isym > 3 || icase > MAXINT)
        break;
      bad_dip = 2;
      if (iv1 >= glob.nvib)
        break;
      if (iv2 >= glob.nvib)
        break;
      if (glob.oblate == FALSE)
        isym = revsym[isym];
      if (isym == 0)
      {
        ld = 0;
        zfac = 0.009274;
        ioff = 8;
      }
      else
      {
        ld = 1;
        iflg |= MELEC;
        ioff = 0;
      }
      imag = 0;
      kd = isoddk[isym]; /* ld = 0,1; kd = 0, 1, 1, 0 */
      good_case = 4;
      if (icase > 0)
      {
        switch (icase)
        {
        case 1: /* ld = 2,2; kd = 0, 1, 1, 2 */
          ld = 2;
          if (isym == 3)
            kd = 2;
          break;
        case 2: /* ld = 2,2; kd = 2, 1, 1, 2 */
        case 11:
        case 12:
          ld = 2;
          if (kd == 0)
            kd = 2;
          if (icase == 11)
          {
            iflg |= MINOQ;
            ++imag;
            good_case = isym; /* commutator */
          }
          else if (icase == 12)
          {
            ifc = -1;
            ++imag;
            isym = 3 - isym;
          }
          break;
        case 5: /* ld = 1,2; kd = 0, 1,  1, 0 */
          if (isym == 0)
          {
            iflg |= MLZ;
            ld = 1;
            ++imag;
            if (lv1 != lv2 || lv1 <= 0 || si1 > 0)
              good_case = 0;
          }
          else
          {
            iflg |= MINOQ;
            ++imag;
            good_case = isym; /* commutator */
          }
          break;
        case 6: /* ld = 0,3; kd = 2, 3, 3, 2 */
          if (isym != 0)
            ld = 3;
          kd += 2;
          break;
        case 7:
          break;
        case 8:
          ifc = -1 - ifc;
          ++imag;
          isym = 3 - isym;
          break;
        case 10:
          iflg |= MINOQ;
          good_case = isym; /* double commutator */
          break;
        }
      }
      bad_dip = 3;
      if (good_case == 0)
        break;
      imag += ld;
      nsym = setgsym((int)pvib1->gsym);
      bad_dip = 4;
      if (testwt(pvib1, pvib2, isym, 0))
        break;
      ldel = lv1 - lv2;
      if (ldel != 0)
      {
        bad_dip = 3;
        if (icase == 5 && isym == 0)
          break;
        if (ldel < 0)
          ldel = -ldel;
        if (ivdif < 0) /* DELTA K * DELTA L < 0 */
          ldel = -ldel;
      }
      else if (lv1 != 0)
      {
        if (glob.newlz)
        {
          if (lv1 < 0)
          {
            iflg |= MLZ;
            isym = 3 - isym;
            ++imag;
          }
        }
        else
        {
          if (lv1 < 0)
            break;
        }
      }
      if (ODD(imag))
        iflg |= MODD;
      /* check for imaginary operators */
      if (iv1 == iv2 && ldel == 0)
      {
        if (ODD(imag))
        {
          if (isym == 0)
            break;
          nimag[isym + ioff] += 1;
          iflg |= MIMAG;
        }
        else if (isym != 0)
        {
          nimag[isym + ioff + 4] += 1;
          iflg |= MIMAG;
        }
      }
      if (glob.nitot >= 3)
      {
        bad_dip = 5;
        if (MOD(kd - ldel, glob.nitot) != 0)
          break; /* require alpha == 0 */
        if (si1 <= itsym)
          iflg |= MIDEN; /* identity operator */
      }
      lv = 1;
      if (si1 > 0)
      {
        /* check spin properties */
        bad_dip = 6;
        if (ld == 0 && isym != 0)
          break;
        iiv1 = pvib1->spt;
        iiv2 = pvib2->spt;
        if (checksp(TRUE, si1, 0, iiv1, iiv2, &zfac) != 0)
          break;
        lv = TEST(iflg, MELEC) ? 2 : 0; /*  simple electric dipole */
      }
      if (ifc < 0)
        iflg |= MFCM;
      if (glob.nofc && ifc != 0)
      {
        if (ifc < 0)
          ifc = -1 - ifc;
        kavg = ifc;
        ifc = 0;
        bad_dip = 7;
        if (kavg > pvib1->knmax && kavg > pvib2->knmax)
          break;
      }
      if (ld == kd && ld < lv)
        ld -= (ld - lv) & 0x7fe; /* subtract even number */
      dipinfo[i].fac = zfac;
      dipinfo[i].flg = (short)iflg;
      dipinfo[i].kd = (signed char)kd;
      dipinfo[i].ld = (signed char)ld;
      dipinfo[i].ldel = (signed char)ldel;
      dipinfo[i].fc = (signed char)ifc;
      dipinfo[i].kavg = (char)kavg;
      bad_dip = 0;
    } while (FALSE);
    if (bad_dip != 0)
    {
      if (bad_dip > 0)
      {
        putbcd(sbcd, NSBCD, &idip[ibcd]);
        if (bad_dip * sizeof(char *) >= sizeof(bad_msg))
          bad_dip = 0;
        fprintf(lu,
                " WARNING: dipole %3d %s has no matrix elements, %s",
                (i + 1), sbcd, bad_msg[bad_dip]);
      }
      memcpy(&dipinfo[i], &dipinfo0, sizeof(SDIP));
      iv2 = iv1 = glob.vibfac;
      isym = 0;
    }
    isimag[i] = 0;
    ii = iv2 * ifac + iv1;
    ijv = ((unsigned int)ii << 2) + (unsigned int)isym;
    if (ndecv < 3)
    {
      idip[ibcd + 5] = idip[ibcd + ndecv + 2];
      idip[ibcd + 4] = idip[ibcd + ndecv + 1];
      idip[ibcd + 3] = (bcd_t)0;
    }
    else
    {
      idip[ibcd + 4] &= 0xf0;
      idip[ibcd + 3] = (bcd_t)((ijv >> 16) & 0xff);
    }
    idip[ibcd + 2] = (bcd_t)(ijv >> 8);
    idip[ibcd + 1] = (bcd_t)ijv;
    if (si1 == 0)
      idip[ibcd + 4] &= 0x0f;
  }
  kmin = 0;
  nmin = ndip;
  imag = 0;
  for (kd = 0; kd < 8; ++kd)
  {
    ld = kd ^ glob.stdphase;
    iv1 = 0;
    iv2 = 0;
    nn = 0;
    for (isym = 1; isym <= 3; ++isym)
    {
      i = ipwr2[isym];
      if (TEST(glob.phasemask, i) && TEST(kd, i))
        break;
      ii = isym;
      if (TEST(ld, i))
        ii += 4;
      nn += nimag[ii + 8];
      iv1 += nimag[ii ^ 4];
      iv2 += nimag[ii];
    }
    if (isym <= 3)
      continue;
    if (iv1 < iv2)
    {
      nn += iv1;
      ii = 1;
    }
    else
    {
      nn += iv2;
      ii = 0;
    }
    if (nn < nmin)
    {
      nmin = nn;
      kmin = kd;
      imag = ii;
      if (nmin == 0)
        break;
    }
  }
  glob.stdphase = (kmin ^ glob.stdphase) & 7;
  if (kmin != 0)
  { /* update phases */
    fprintf(lu,
            "NON-STANDARD PHASE CONVENTION IN USE, %2d\n", glob.stdphase);
    for (isym = 1; isym <= 3; ++isym)
    {
      ixphase[isym] = (ixphase[isym] + kmin) & 1;
      kmin = kmin >> 1;
    }
  }
  if (imag != 0)
  { /* fix flags */
    for (i = 0; i < ndip; ++i)
    {
      iflg = dipinfo[i].flg;
      if ((iflg & MELEC) == 0)
        continue;
      iflg |= MDIPI;
      dipinfo[i].flg = (short)iflg;
      isimag[i] = 1;
    }
  }
  if (nmin > 0)
  {
    for (i = 0; i < ndip; ++i)
    {
      iflg = dipinfo[i].flg;
      if ((iflg & MIMAG) == 0)
        continue;
      isym = (int)(idip[ibcd + 1] & 3);
      if (isym == 0)
        continue;
      ii = ixphase[isym] + isym;
      if (kmin != 0 && TEST(iflg, MELEC))
        ++ii;
      if (TEST(iflg, MODD))
        ++i;
      if (EVEN(ii))
        continue;
      /* set dipole to zero */
      putbcd(sbcd, NSBCD, &idip[ibcd]);
      fprintf(lu,
              " WARNING: dipole %3d %s has bad phase and is ignored\n",
              (i + 1), sbcd);
      dipinfo[i].fac = 0.;
      dipinfo[i].flg = (short)0;
      iv1 = glob.vibfac;
      isym = 0;
      ii = iv1 * ifac + iv1;
      ijv = ((unsigned int)ii << 2) + (unsigned int)isym;
      if (ndecv < 3)
      {
        idip[ibcd + 3] = (bcd_t)0;
      }
      else
      {
        idip[ibcd + 4] &= 0xf0;
        idip[ibcd + 3] = (bcd_t)((ijv >> 16) & 0xff);
      }
      idip[ibcd + 2] = (bcd_t)(ijv >> 8);
      idip[ibcd + 1] = (bcd_t)ijv;
      if (--nmin == 0)
        break;
    }
  }
  return 0;
} /* setint */

/**
 * @brief Calculates intensity matrix elements for transitions between states in different blocks.
 * @param iblk_upper_state Upper state F-block index.
 * @param isiz_upper_dim Dimension of upper state F-block.
 * @param jblk_lower_state Lower state F-block index.
 * @param jsiz_lower_dim Dimension of lower state F-block.
 * @param num_dipole_params Number of dipole moment parameters.
 * @param bcd_dipole_ids Array of dipole moment identifiers in BCD format.
 * @param dipole_values Array of dipole moment values.
 * @param intensity_matrix_out Output intensity matrix S_ij = <bra_i|mu|ket_j>.
 * @return int (2*F_upper + 1) if transitions are allowed and calculated, 0 otherwise.
 *         Effectively, the degeneracy of the upper F state.
 */
int intens(iblk, isiz, jblk, jsiz, ndip, idip, dip, s)
    const int iblk,
    isiz, jblk, jsiz, ndip;
const bcd_t *idip;
const double *dip;
double *s;
{
  /*  function to calculate intensities */
  /*     on entry: */
  /*         IBLK,JBLK= block numbers */
  /*         ISIZ,JSIZ= block sizes */
  /*         IDIP= dipole id codes , IDIP(0)= number of dipoles */
  /*         DIP= dipole values */
  /*     on return: */
  /*         S= transition dipole matrix in original basis */
  /*         INTENS returns true is S is not 0 */
  static int idipoff;
  SDIP *pdip;
  double dd;
  long ndms;
  unsigned int ijv;
  int ncos, ifup, isym, i, k, lv, n, icase, ibase, mkd, ibcd, nbcd, isunit;
  int jbase, kbgni, kbgnj, nblki, nblkj, ixtmp, iret, kd, ld, npair, ldel;
  int ix, jx, nx, si1, si2, iff, jff, nni, nnj, iwt[5], jwt[5], ndecv;
  int ivv, jvv, ixx, jxx, ifc, ksym, kl, ioff, joff, alpha, nitot, dipoff;
  bcd_t bijv1, bijv2, bijv3;

  if (ndmx <= 0)
  {
    puts("working vectors not allocated");
    exit(EXIT_FAILURE);
  }
  dipoff = idipoff;
  if (dipoff >= nddip)
    dipoff = 0;
  idipoff = dipoff + ndip;
  iret = 0;
  cgetv[0].cblk = 0;
  cgetv[1].cblk = 0;
  /*     get quantum information */
  nblki = getqq(iblk, &iff, iwt, ibkptr, ikmin, ivs);
  ioff = glob.maxblk;
  nblkj = getqq(jblk, &jff, jwt, &ibkptr[ioff], &ikmin[ioff], &ivs[ioff]);
  ixx = iff - jff;
  if (ixx < -2 || ixx > 2)
    return 0;
  ixx = iff + jff;
  if (ixx < 2 || ODD(ixx))
    return 0;
  if (iwt[3] != jwt[3])
    return 0;
  if (checkwt(iwt, jwt) != 0)
    return 0;
  nitot = glob.nitot;
  alpha = 0;
  bijv3 = (bcd_t)0;
  /* clear dipole matrix */
  dclr(isiz, jsiz, s, 1);
  /* loop over sub-blocks */
  ndms = isiz;
  ndecv = glob.vibdec;
  nbcd = (int)idip[0] & 0x7f;
  for (ixx = 0; ixx < nblki; ++ixx)
  {
    ibase = ibkptr[ixx];
    n = ibkptr[ixx + 1] - ibase;
    kbgni = ikmin[ixx];
    ixtmp = ivs[ixx];
    getqs(ixtmp, iff, n, kbgni, ixcom, iscom, &ivv);
    nni = iscom[0];
    for (jxx = 0; jxx < nblkj; ++jxx)
    {
      joff = jxx + ioff;
      jbase = ibkptr[joff];
      n = ibkptr[joff + 1] - jbase;
      kbgnj = ikmin[joff];
      ixtmp = ivs[joff];
      getqs(ixtmp, jff, n, kbgnj, jxcom, jscom, &jvv);
      nnj = jscom[0];
      nx = ((nni > nnj) ? nni : nnj) >> 1;
      i = (ivv < jvv) ? ivv : jvv;
      ijv = (unsigned int)(ivv + jvv + i * glob.vibfac);
      ijv = (ijv << 2) + blksym(ixcom, jxcom);
      bijv1 = (bcd_t)ijv;
      bijv2 = (bcd_t)(ijv >> 8);
      if (ndecv == 3)
        bijv3 = (bcd_t)(ijv >> 16);
      for (i = 0, ibcd = 0; i < ndip; ++i, ibcd += nbcd)
      {
        /*  set dipole for transition */
        if (bijv1 != idip[ibcd + 1])
          continue;
        if (bijv2 != idip[ibcd + 2])
          continue;
        if (bijv3 != idip[ibcd + 3])
          continue;
        pdip = &dipinfo[dipoff + i];
        kl = pdip->flg;
        if (TEST(kl, MINOQ) && nni == nnj)
          continue;
        isym = (int)(bijv1 & 3);
        si1 = (int)(idip[ibcd + 4] >> 4) & 0x0f;
        si2 = 0;
        icase = bcd2i(idip[ibcd + 5]);
        kd = (int)pdip->kd;
        ifc = (int)pdip->fc;
        ldel = (int)pdip->ldel;
        lv = 1;
        ld = (int)pdip->ld;
        if (si1 > 0)
          lv = TEST(kl, MELEC) ? 2 : 0;
        ksym = 0;
        if (nitot >= 3 && kd > 0)
        {
          alpha = nitot;
          ksym = 1;
        }
        do
        {
          if (ksym == 0)
          {
            alpha = 0;
          }
          /*  find matrix elements */
          mkd = getmask(ixcom, jxcom, kd, ldel, kl, alpha);
          if (mkd == 0)
            continue;
          npair = getll(2, ld, lv, kd, si1, si2, lscom, iscom, jscom);
          if (npair < 0)
            continue;
          ifup = TEST(kl, MDIPI) ? -1 : 1;
          isunit = (int)pdip->kavg;
          if (isunit > nx)
            continue;
          ncos = dircos(ixcom, jxcom, ld, kd, ndmx, wk, idx, jdx,
                        ifup, kl, mkd, &isunit);
          if (ncos == 0)
            continue;
          if (ncos < 0)
          {
            puts("DIRCOS WORKING VECTOR TOO SHORT IN INTENS");
            exit(EXIT_FAILURE);
          }
          dd = dip[i] * pdip->fac;
          /* correct for special cases */
          switch (icase)
          {
          case 3:
            if (nni == nnj)
            {
              k = nx * (nx + 1);
            }
            else
            {
              k = nx * nx;
            }
            dd *= k;
            break;
          case 4:
            symksq(1, kbgni, kbgnj, ncos, wk, idx, jdx);
            break;
          case 5:
          case 11:
            if (ifup > 0)
              break;
            if (nni < nnj)
              dd = -dd;
            dd *= nx;
            break;
          case 9:
            if (nni == nnj)
            {
              k = nx * (nx + 1);
              dd *= k;
            }
            else
            {
              k = nx * nx;
              dd *= k + 1;
            }
            dd *= k;
            break;
          case 10:
            dd *= nx * nx;
            break;
          case 7:
          case 8:
          case 12:
            specfc(ifc, ivv, jvv, kd, kbgni, kbgnj, ncos, wk, idx, jdx);
            break;
          }
          /*  correct for reduced matrix of N */
          if (ld != lv)
            dd *= rmatrx(ld, lv, ixcom, jxcom);
          /*  couple dipoles through the spins */
          tensor(&dd, iscom, jscom, lscom, ismap, npair, alpha);
          for (k = 0; k < ncos; ++k)
          {
            ix = idx[k] + ibase;
            jx = jdx[k] + jbase;
            s[ix + jx * ndms] += dd * wk[k];
          }
          iret = iff + 1;
        } while (--ksym >= 0);
      }
    }
  }
  return iret;
} /* intens */

/**
 * @brief Reads and sets up global options for calculations from an input file.
 *        This includes vibrational state definitions, spin configurations, K ranges,
 *        statistical weights, symmetry options, and diagonalization preferences.
 * @param input_file_unit File pointer for the input option file.
 * @param num_qn_formats_out Output: Pointer to store the number of quantum number formats determined.
 * @param molecule_type_out Output: Pointer to store molecule type (2 for linear, 3 for non-linear). Renamed itd.
 * @param num_bcd_digits_param_id_out Output: Pointer to store number of BCD digits per parameter ID. Renamed nbcd.
 * @param parameter_name_file_path_out Output: Character array to store the name of the parameter name file (e.g., "sping.nam"). Renamed namfil.
 * @return int Number of option cards read from the input file, or -1 if end-of-file is reached prematurely.
 */
int setopt(luin, nfmt, itd, nbcd, namfil)
FILE *luin;
int *nfmt, *itd, *nbcd;
char *namfil;
{
  /* this subroutine reads in quantum quantities */
  /*     on entry: */
  /*         LUIN is the fortran unit number to read from */
  /*         NFMT is the maximum number of quanta */
  /*     on return: */
  /*         NFMT= number of quantum number formats for catalog */
  /*         ITD= 2 for linear, 3 for non-linear */
  /*         NAMFIL contains a file name on which to find parameter names */
  /*         SETOPT = number of cards read, -1 for end-of-file */
#define NVEC 11
#define NDECSPIN 11
#define MAXIAX 11
#define NCARD 130
  static int isymv[] = {2, 2, 2, 2, 2, 2, 3, 4, 4, 5, 6, 6};
  static int iaxv[] = {2, 1, 2, 3, 4, 5, 2, 4, 5, 2, 4, 5};
  SVIB *pvinfo, *pvinfom;
  short *pshort;
  double rvec[NVEC], rvec0[NVEC], vsym;
  size_t nl;
  int ivib, lopt, i, ii, j, k, nvibm, isym, iend, ivsym;
  int iwtpl, iwtmn, nt, ncards, knnmin, knnmax, si2, iax;
  int spinsgn, ewt0, ewt1, nspinqn;
  int ivwt[3];
  bcd_t bcdspin[NDECSPIN];
  char card[NCARD], ctyp;

  cgetv[0].cblk = cgetv[1].cblk = 0; /* reset store for getqn */
  cjjini();
  glob.lsym = TRUE;
  glob.esym = TRUE;
  glob.esymdec = 100;
  glob.maxwt = 0;
  glob.stdphase = 0;
  glob.newlz = FALSE;
  glob.nofc = FALSE;
  glob.g12 = 0;
  rvec0[0] = 1;
  rvec0[1] = 0;
  rvec0[2] = (double)MAXN_DIRCOS;
  rvec0[3] = 0;
  rvec0[4] = 1;
  rvec0[5] = 1;
  rvec0[6] = 1;
  rvec0[7] = 0;
  rvec0[8] = 99;
  rvec0[9] = 0;
  rvec0[10] = 0;
  ncards = nvibm = nt = si2 = nspinqn = 0;
  *itd = 2;
  bcdspin[0] = (bcd_t)NDECSPIN;
  /*     read option cards */
  do
  {
    if (fgetstr(card, NCARD, luin) <= 0)
      return -1;
    ++ncards;
    ctyp = C0;
    if (isalpha((int)card[0]))
      ctyp = card[0];
    else if (isalpha((int)card[1]))
      ctyp = card[1];
    iend = getbcd(card, bcdspin, NDECSPIN);
    spinsgn = NEGBCD(bcdspin[0]);
    dcopy(NVEC, rvec0, 1, rvec, 1);
    rvec[7] = 0.;
    pcard(&card[iend], rvec, NVEC, NULL);
    lopt = (int)rvec[0];
    knnmin = (int)rvec[1];
    knnmax = (int)rvec[2];
    iax = (int)rvec[4];
    iwtpl = (int)rvec[5];
    iwtmn = (int)rvec[6];
    vsym = rvec[7];
    ewt0 = (int)rvec[8];
    ewt1 = 0;
    if (vsym < -0.5)
    {
      ivsym = -1;
      vsym = 0.;
    }
    else if (vsym > 0.5 && vsym < 2.e+15)
    {
      ivsym = 1;
    }
    else
    {
      ivsym = 0;
      vsym = 0.;
    }
    ivib = lopt;
    if (ivib < 0)
      ivib = -ivib;
    if (knnmax < 0 || knnmax > MAXN_DIRCOS)
      knnmax = MAXN_DIRCOS;
    if (knnmin > knnmax)
      knnmax = knnmin;
    if (knnmin != knnmax)
      *itd = 3;
    if (ncards == 1)
    {
      dcopy(NVEC, rvec, 1, rvec0, 1);
      strcpy(namfil, "sping.nam");
      if (ctyp != C0)
        namfil[4] = ctyp;
      glob.ixz = (int)rvec[3];
      glob.idiag = (int)rvec[9];
      k = (int)rvec[10];
      glob.stdphase = MOD(k, 10);
      k = k / 10;
      if (ODD(k))
        glob.newlz = TRUE;
      if (ODD2(k))
      {
        glob.nofc = TRUE;
      }
      else if ((k & 4) != 0)
      {
        glob.g12 = 1;
      }
      if (ivib > 80)
      {
        if (ivib > MAXVIB)
          ivib = MAXVIB;
        if (ivib > 80 && sizeof(int) == 2)
          ivib = 80;
      }
      else if (ivib <= 0)
      {
        ivib = 1;
      }
      if (glob.nvib > 1)
      {
        free(vinfo);
        vinfo = NULL;
      }
      vinfo1.spt = sptzero;
      if (ivib > 1)
      {
        nvibm = ivib - 1;
        nl = (size_t)ivib * sizeof(SVIB);
        vinfo = (SVIB *)mallocq(nl);
        memcpy(vinfo, &vinfo1, sizeof(SVIB));
      }
      else
      {
        vinfo = &vinfo1;
      }
      for (k = 0; k <= 4; ++k)
        vinfo->wt[k] = 0;
      vinfo->ewt[0] = vinfo->ewt[1] = 0;
      glob.nvib = ivib;
      if (ivib <= 9)
      {
        glob.vibfac = 9;
        glob.vibdec = 1;
        k = 4; /* 16 */
      }
      else if (ivib <= 99)
      {
        glob.vibfac = 99;
        glob.vibdec = 2;
        k = 7; /* 128 */
      }
      else
      {
        glob.vibfac = 999;
        glob.vibdec = 3;
        k = 10; /* 1024 */
      }
      glob.msshft = (unsigned int)(k + 2);
      /* msshft = shifts needed for vib and symmetry */
      glob.msmask = (int)(1 << (unsigned int)k) - 1;
      ivib = 0;
      glob.oblate = (lopt < 0);
      glob.nqnn = 3;
      if (spinsgn != 0)
        glob.nqnn = 1;
      if (ewt0 < 0)
        glob.esymdec = 1000;
    }
    else if (ivib >= glob.nvib)
    {
      /* if IVIB out of range read another card */
      continue;
    }
    pvinfo = &vinfo[ivib];
    k = iax;
    if (iax < 0)
      iax = -iax;
    if (iax > MAXIAX)
      iax = 0;
    if (iax <= 3 && glob.oblate)
      iax = revsym[iax];
    isym = isymv[iax];
    iax = iaxv[iax];
    pvinfo->gsym = (short)(isym << 1);
    if (isym >= 3)
    {
      j = 2 - (isym & 1);
      if (glob.maxwt < j)
        glob.maxwt = j;
    }
    /*  set up minimum K values */
    pvinfo->knmax = (short)knnmax;
    pvinfo->knmin[0] = 0;
    pvinfo->knmin[1] = 1;
    pvinfo->knmin[2] = 1;
    pvinfo->knmin[3] = 2;
    if (k < 0)
    { /* iax was negative */
      if (isym >= 3)
      {
        pvinfo->knmin[0] = -2;
        pvinfo->knmin[3] = -2;
      }
      pvinfo->gsym |= 1;
    }
    if (knnmax != 0 && glob.nqnn == 1)
      glob.nqnn = 2;
    /* set ewt */
    if (ewt0 < 0)
      ewt0 = -ewt0;
    if (ewt0 >= glob.esymdec)
    {
      ewt1 = ewt0 / glob.esymdec;
      ewt0 -= glob.esymdec * ewt1;
    }
    if (isym < 3)
    {
      ewt0 = glob.esymdec - 1;
      glob.esym = FALSE;
    }
    else if (isym == 4 && iax == 5)
    { /* 4-fold B entry has no E */
      ewt0 = -1;
    }
    if (isym > 3 && pvinfo->wt[4] < 4)
    {
      pvinfo->ewt[0] = pvinfo->ewt[1] = (short)0;
      pvinfo->wt[4] |= 4;
    }
    if (isym == 6)
    { /* six-fold, iax = 4,5 */
      pvinfo->ewt[iax - 4] = (short)ewt0;
    }
    else if (ewt0 >= 0)
    {
      pvinfo->ewt[0] = pvinfo->ewt[1] = (short)ewt0;
    }
    pvinfo->lvupper = 0;
    if (ewt1 != 0)
    {
      if (ewt1 > 10)
      {
        pvinfo->lvupper = 2;
        ewt1 -= 10;
      }
      k = isym;
      if (k < 3)
        k = 3;
      if (ewt1 > k)
        ewt1 = MOD(ewt1 - 1, k) + 1;
      pvinfo->lvupper += (short)(ewt1 << 2);
      pvinfo->knmin[0] = 0;
      pvinfo->knmin[3] = 0; /* 0,1,1,0 */
    }
    pvinfo->lvqn = (short)ewt1;
    for (i = 0; i <= 3; ++i)
    {
      k = pvinfo->knmin[i];
      if (k < 0 && knnmin == 0)
        continue;
      while (k < knnmin)
        k += 2;
      pvinfo->knmin[i] = (short)k;
    }
    k = getsp(bcdspin, pvinfo);
    if (nspinqn < k)
      nspinqn = k;
    pvinfo->nqn = (short)k;
    if (ncards == 1)
    { /* fill in defaults with wt = 0 */
      pvinfom = vinfo;
      for (i = 1; i <= nvibm; ++i)
      {
        ++pvinfom;
        memcpy(pvinfom, vinfo, sizeof(SVIB));
      }
      ivib = -glob.nvib;
    }
    setwt(pvinfo, ivib, iax, iwtpl, iwtmn, vsym); /* set up weights */
  } while (ivsym < 0); /*   end of reading */
  setsp();
  /* set up format */
  if (glob.nqnn == 3)
    *itd = 3;
  glob.nqn = glob.nqnn;
  glob.vibfmt = (nvibm != 0);
  glob.iqfmt0 = 0;
  if (glob.vibfmt)
  {
    glob.nqn += 1;
    glob.iqfmt0 = 1000;
  }
  nt = glob.nqn;
  k = nt + nspinqn;
  glob.nqn = k;
  if (k > (*nfmt))
    k = nt + 2;
  glob.maxqn = k;
  glob.nqn0 = nt;
  glob.iqfmt0 += 100 * nt + k;
  /*     make sure +/- l values come in pairs */
  pvinfo = pvinfom = vinfo;
  k = 0;
  for (i = 0; i <= nvibm; ++i)
  {
    if (i > 0)
    {
      pvinfom = pvinfo;
      ++pvinfo;
    }
    pvinfo->nqn += (short)nt;
    pvinfo->wt[4] &= 3;
    j = pvinfo->wt[4];
    if (EVEN(pvinfo->gsym))
    {
      pvinfo->wt[4] = 0; /* no Itot */
    }
    else if (pvinfo->wt[0] == pvinfo->wt[j] &&
             pvinfo->wt[1] == pvinfo->wt[j ^ 1])
    {
      pvinfo->wt[4] = (short)((pvinfo->gsym > 6) ? -1 : 0);
    }
    j = pvinfo->lvqn;
    if (k != 0)
    { /* previous l was non-zero and first of pair */
      isym = (int)((unsigned int)pvinfo->gsym >> 1);
      if (isym >= 3)
      {
        if ((k + k) == isym && j != k)
        {        /* B symmetry */
          k = j; /* single instance */
          pvinfom->lvqn = 0;
          pshort = pvinfom->wt;
          ii = pshort[0];
          pshort[0] = pshort[2];
          pshort[2] = (short)ii;
          ii = pshort[1];
          pshort[1] = pshort[3];
          pshort[3] = (short)ii;
          pshort = pvinfom->ewt;
          ii = pshort[0];
          pshort[0] = pshort[1];
          pshort[1] = (short)ii;
          continue;
        }
        else
        {
          if (MOD(j + k, isym) != 0)
            break;
        }
      }
      if (k > j)
      {
        if (j == 0)
          break;
        pvinfom->lvqn = (short)(-j);
        if (pvinfo->knmin[0] == 0)
          pvinfo->knmin[0] = 2;
        if (pvinfo->knmin[3] == 0)
          pvinfo->knmin[3] = 2;
      }
      else
      {
        pvinfo->lvqn = (short)(-k);
        if (pvinfom->knmin[0] == 0)
          pvinfom->knmin[0] = 2;
        if (pvinfom->knmin[3] == 0)
          pvinfom->knmin[3] = 2;
      }
      for (j = 0; j <= 4; ++j)
      {
        if (pvinfo->wt[j] != pvinfom->wt[j])
          break;
      }
      if (j <= 4)
        break;
      if (pvinfo->gsym != pvinfom->gsym)
        break;
      if (pvinfo->ewt[0] != pvinfom->ewt[0])
        break;
      if (pvinfo->ewt[1] != pvinfom->ewt[1])
        break;
      pvinfo->lvupper += 1;
      j = 0;
    }
    k = j;
  }
  if (k != 0)
  {
    puts("L DOUBLETS SHOULD BE IN PAIRS");
    exit(EXIT_FAILURE);
  }
  k = 0;
  pvinfo = vinfo;
  for (ivib = 0; ivib <= nvibm; ++ivib)
  {
    nt = pvinfo->nspstat;
    ii = setgsym((int)pvinfo->gsym);
    for (isym = 0; isym < 4; ++isym)
    {
      if (getwt(pvinfo, isym, 0, ivwt) <= 0)
        continue;
      for (ii = 1; ii <= nt; ++ii)
      {
        if (getwt(pvinfo, isym, ii, ivwt) > 0)
          ++k;
      }
    }
    ++pvinfo;
  }
  if (glob.nbkpj > 0)
  {
    free(moldv);
    moldv = NULL;
    free(blkptr);
    blkptr = NULL;
  }
  glob.nbkpj = k;
  nl = (size_t)k * sizeof(int);
  moldv = (int *)mallocq(nl);
  moldv[0] = 0;
  nl += sizeof(int);
  blkptr = (int *)mallocq(nl);
  blkptr[0] = 0;
  k = 0;
  pvinfo = vinfo;
  for (ivib = 0; ivib <= nvibm; ++ivib)
  {
    nt = pvinfo->nspstat;
    ii = setgsym((int)pvinfo->gsym);
    for (isym = 0; isym < 4; ++isym)
    {
      if (getwt(pvinfo, isym, 0, ivwt) <= 0)
        continue;
      for (ii = 1; ii <= nt; ++ii)
      {
        if (getwt(pvinfo, isym, ii, ivwt) > 0)
        {
          blkptr[k] = k;
          moldv[k] = isym + (ivib << 2) + (ii << glob.msshft);
          ++k;
        }
      }
    }
    ++pvinfo;
  }
  blkptr[k] = k;
  pvinfo = NULL;
  *nbcd = (NDECPAR + 1) + glob.vibdec;
  if (sizeof(int) == 2 && glob.vibdec > 2)
  {
    puts("too many vibrations for 16-bit integer computer");
    ncards = 0;
  }
  *nfmt = 1;
  if (nspinqn != 0 && nvibm != 0)
    *nfmt = setfmt(&k, -1);
  return ncards;
} /* setopt */

/**
 * @brief Determines the quantum number output format code (QNFMT) for vibrational states.
 *        If nfmt_or_vib_limit > 0, it fills iqnfmt_out array with format codes for each vibrational state.
 *        If nfmt_or_vib_limit = -1, it checks if all states up to glob.nvib have the same format code
 *        and returns 1 if they do, or glob.nvib if they differ (this signals multiple formats needed).
 * @param iqnfmt_out Output: If nfmt_or_vib_limit > 0, this array is filled with QNFMT codes for each vibrational state.
 *                   If nfmt_or_vib_limit = -1, this parameter is effectively a dummy (value passed as &k from setopt).
 * @param nfmt_or_vib_limit If > 0, it's the limit of vibrational states to process for filling iqnfmt_out.
 *                          If -1, it's a flag to check format consistency across all vibrational states.
 * @return int If nfmt_or_vib_limit > 0, returns glob.maxqn (total QN fields to be output).
 *             If nfmt_or_vib_limit = -1, returns 1 if all formats are same, or a value > 1 (e.g. num_vib_states_to_check) if they differ.
 */
int setfmt(iqnfmt, nfmt)
int *iqnfmt, nfmt;
{
  SVIB *pvinfo;
  short *jjs;
  int iqfmtq, iqfmt0, i, j, k, iv, si2, ff, nset, nspinv, nvib;
  iqfmt0 = 0;
  nvib = glob.nvib;
  if (nfmt > 0 && nfmt < nvib)
    nvib = nfmt;
  pvinfo = vinfo;
  for (iv = 0; iv < nvib; ++iv)
  {
    iqfmtq = glob.iqfmt0;
    i = setgsym((int)pvinfo->gsym);
    jjs = pvinfo->spt;
    nset = jjs[0];
    nspinv = nset - 1;
    if (nspinv > nspin)
      nspinv = nspin;
    jjs += nset;
    ff = 200 - jjs[0];
    if ((int)pvinfo->nqn > glob.maxqn)
    {
      iqfmtq += 10 * (ff & 1) + 4000;
    }
    else
    {
      si2 = 0;
      k = glob.maxqn - glob.nqn0 - 1;
      for (i = 0; i < k; ++i)
      {
        if (i == itsym && itsym < itptr)
        {
          i = itptr;
          j = 0;
          si2 = si2 << 1;
        }
        if (i < nspinv)
        {
          j = jjs[i + 1];
          if (i < itsym)
            j += ff;
        }
        else
        {
          j = ff;
        }
        si2 = (si2 + (j & 1)) << 1;
      }
      si2 += (ff & 1);
      if (si2 > 9)
        si2 &= 7;
      if (itptr < nspin)
        iqfmtq += 2000;
      if (itsym < itptr)
        iqfmtq += 4000;
      iqfmtq += 10 * si2;
    }
    if (iv == 0)
      iqfmt0 = iqfmtq;
    if (nfmt < 0)
    {
      if (iqfmtq != iqfmt0)
        return nvib;
    }
    else
    {
      iqnfmt[iv] = iqfmtq;
    }
    ++pvinfo;
  }
  if (nfmt <= 0)
    return 1;
  return glob.maxqn;
} /* setfmt */

/**
 * @brief Initializes global parameters and allocates/frees memory related to Hamiltonian blocks
 *        and parameter processing. Called from main program to set up before calculations or
 *        to clean up.
 * @param list_unit File pointer for listing output. Renamed lu.
 * @param num_total_params Total number of parameters. Renamed npar.
 * @param bcd_parameter_ids Array of BCD parameter identifiers. Renamed idpar.
 * @param parameter_values Array of parameter values. Renamed par.
 * @param max_f_qn_in_block Pointer to store maximum F quantum number used in any block. Renamed nbkpf.
 * @param max_hamiltonian_dim Pointer to store maximum Hamiltonian dimension. Renamed negy.
 * @return int Returns a value related to vibrational state packing (vibfac^2 related), or 0 if npar is 0.
 *         The exact meaning of the return value might be specific to the calling context.
 */
int setblk(lu, npar, idpar, par, nbkpf, negy)
FILE *lu;
const int npar;
int *nbkpf, *negy;
bcd_t *idpar;
const double *par;
{
  /* this subroutine initializes  quantum quantities */
  /*     on entry: */
  /*         LU= logical unit for listing */
  /*         IDPAR= parameter id */
  /*         NBKPF= maximum value of F */
  /*         NEGY = max dimensioned energy in main program */
  /*     on return: */
  /*         NBKPF= blocks per F */
  /*         NEGY= possibly downsized NEGY to account for local storage */
  static char csym[] = {'x', 'c', 'b', 'a'};
  SVIB *pvinfo;
  SPAR *spar_now;
  short *jjs;
  double ai;
  size_t nl;
  int nodd, isym, jsym, ixzt, i, k, kl, lt, n, kdiag, nvibm, neven, nsize, alpha;
  int kd, ld, ii, jj, im, kk, ni, nj, iv, jv, ivmin, ix, lv, lx, npair, ldel;
  int knnmin, ntstat, si1, si2, iff, ikq, jjt, ins, neuler, sznz, iff0;
  int iwt, nsym, gsym, maxblk, maxsblk, itmp, nd, ifc, nf, kavg;
  unsigned int ivsym;
  int ivwt[3], jvwt[3];
  char ctmp;

  cgetv[0].cblk = cgetv[1].cblk = 0; /* reset store for getqn */
  kdiag = glob.idiag;
  if (npar <= 0)
  {
    if (glob.nvib > 1)
      free(vinfo);
    vinfo = &vinfo1;
    glob.nvib = 0;
    if (glob.nbkpj > 0)
    {
      free(moldv);
      moldv = &zmoldv;
      free(blkptr);
      blkptr = &zblkptr;
      glob.nbkpj = 0;
    }
    if (glob.maxblk > 0)
    {
      free(ivs);
      ivs = &zivs;
      free(ikmin);
      ikmin = &zikmin;
      free(ibkptr);
      ibkptr = &zibkptr;
      glob.maxblk = 0;
    }
    if (ndmx > 0)
    {
      free(iqnsep);
      iqnsep = &ziqnsep;
      free(idx);
      idx = &zidx;
      free(jdx);
      jdx = &zjdx;
      free(wk);
      wk = &zwk;
      ndmx = 0;
    }
    return 0;
  }
  nvibm = glob.nvib - 1;
  n = (*nbkpf);
  if (n > 0)
  {
    iff0 = n << 1;
    pvinfo = vinfo;
    for (i = 0; i <= nvibm; ++i)
    {
      jjs = pvinfo->spt;
      k = jjs[0];
      /*  select largest N relative to F */
      ni = (iff0 - jjs[k]) >> 1;
      if ((int)pvinfo->knmax > ni)
        pvinfo->knmax = (short)ni;
      ++pvinfo;
    }
    pvinfo = NULL;
  }
  if (pasort(lu, npar, idpar, par) == 0)
    glob.idiag = -1;
  pvinfo = NULL;
  iff0 = 200;
  /*    find interactions */
  ntstat = glob.nbkpj;
  sznz = ndmax = 0;
  if (ntstat <= 0)
    return 0;
  nf = nspin - 1;
  for (ii = 0; ii < ntstat; ++ii)
  {
    im = moldv[ii];
    iff = iff0;
    kk = getqs(im, iff, -1, -1, ixcom, iscom, &iv);
    if (ODD(iscom[nspin]))
    {
      /* special for odd spin multiplicity (make N integer) */
      kk = getqs(im, --iff, 0, -1, ixcom, iscom, &iv);
    }
    isym = ixcom[XSYM];
    ni = ixcom[XNVAL];
    pvinfo = &vinfo[ixcom[XVIB]];
    gsym = pvinfo->gsym;
    nsym = setgsym(gsym);
    getwt(pvinfo, isym, kk, ivwt);
    for (jj = 0; jj < ii; ++jj)
    {
      if (blkptr[ii] == blkptr[jj])
        continue;
      jjt = moldv[jj];
      kk = getqs(jjt, iff, 0, -1, jxcom, jscom, &jv);
      jsym = jxcom[XSYM];
      nj = jxcom[XNVAL];
      nd = ni - nj;
      pvinfo = &vinfo[jxcom[XVIB]];
      if (pvinfo->gsym != (short)gsym)
        continue;
      getwt(pvinfo, jsym, kk, jvwt);
      for (k = glob.maxwt; k >= 0; --k)
      {
        if (ivwt[k] != jvwt[k])
          break;
      }
      if (k >= 0)
        continue;
      /* check for connection restrictions on N */
      ixzt = glob.ixz;
      if (ODD(ixzt) && iscom[nspin] != jscom[nspin])
        continue;
      for (k = 0; k < nf; ++k)
      {
        if (k == itsym)
          k = itptr;
        ixzt = ixzt >> 1;
        /* check for matching spin multiplicity */
        if (ODD(iscom[k] + jscom[k]))
          break;
        /* check for connection restrictions on spin */
        if (ODD(ixzt) && iscom[k] != jscom[k])
          break;
      }
      if (k < nf)
        continue;
      ivmin = (iv < jv) ? iv : jv;
      itmp = iv + jv + ivmin * glob.vibfac;
      ivsym = blksym(ixcom, jxcom) + ((unsigned int)itmp << 2);
      for (spar_now = spar_head[ivmin]; spar_now != NULL;
           spar_now = spar_now->next)
      {
        if (spar_now->ipsym < ivsym)
          break;
        if (spar_now->ipsym > ivsym)
          continue;
        if (TEST(spar_now->flags, MNSQ))
          continue;
        if (sznz < 0)
          sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom);
        kl = idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins, &si1, &si2,
                    &sznz, &ifc, &alpha, &ldel, &kavg);
        if (ODD(neuler))
          continue;
        if (sznz > 0)
        {
          sznz = sznzfix(sznz, ni, nj, ixcom, jxcom, iscom, jscom);
        }
        if (getmask(ixcom, jxcom, kd, ldel, kl, alpha) == 0)
          continue;
        npair = getll(0, ld, lt, kd, si1, si2, lscom, iscom, jscom);
        if (npair < 0)
          continue;
        if (kd > 0)
          glob.idiag = kdiag;
        if (nd > ndmax)
          ndmax = nd;
        kd = blkptr[ii] - blkptr[jj];
        if (kd != 0)
        { /* II and JJ are coupled */
          if (kd > 0)
          { /* rename ptr II to ptr JJ */
            kk = blkptr[ii];
            ix = jj;
          }
          else
          { /* rename ptr JJ to ptr II */
            kk = blkptr[jj];
            ix = ii;
          }
          glob.idiag = kdiag;
          for (k = 0; k <= ii; ++k)
          {
            if (blkptr[k] == kk)
              blkptr[k] = blkptr[ix];
          }
        }
        break;
      } /* loop over parameters */
    } /* jj loop */
  } /* ii loop */
  /*  block connections now found */
  if (glob.idiag < 0)
    glob.idiag = -1;
  switch (glob.idiag)
  {
  case -1:
    fputs("NO DIAGONALIZATION NEEDED\n", lu);
    break;
  default:
    fputs("ENERGY SORT OF WANG SUB-BLOCKS\n", lu);
    break;
  case 1:
    fputs("EIGENVECTOR SORT OF STATES\n", lu);
    break;
  case 2:
    fputs("ENERGY FOLLOWS FIRST ORDER WITHIN WANG SUB-BLOCK\n", lu);
    break;
  case 3:
    fputs("ENERGY FOLLOWS ORDER OF TAU WITHIN V AND SPIN SUB-BLOCKS\n", lu);
    break;
  case 4:
    fputs("ENERGY FOLLOWS K*K WITHIN V AND SPIN SUB-BLOCKS\n", lu);
    break;
  case 5:
    fputs("ENERGY FOLLOWS FIRST ORDER WITHIN V AND SPIN SUB-BLOCKS\n", lu);
  }
  if (glob.stdphase > 0)
    fprintf(lu, "NON-STANDARD PHASE CONVENTION IN USE, %2d\n", glob.stdphase);
  if (glob.newlz)
    fputs("Lz DEFINED EXPLICITLY\n", lu);
  if (glob.nofc)
    fputs("ALTERNATE DEFINITION OF PARAMETER FC FIELD\n", lu);
  if (glob.g12 != 0)
    fputs("G12 group alternating sign of Fourier coefficients with K  \n", lu);
  if (glob.oblate)
  {
    fputs("OBLATE ROTOR\n", lu);
    ctmp = csym[1];
    csym[1] = csym[3];
    csym[3] = ctmp;
  }
  else
  {
    fputs("PROLATE ROTOR\n", lu);
  }
  if (glob.nqnn == 1)
  {
    fputs("LINEAR MOLECULE QUANTA, K SUPPRESSED\n", lu);
  }
  else if (glob.nqnn == 2)
  {
    fputs("SYMMETRIC TOP QUANTA\n", lu);
  }
  fputs("    V KMIN KMAX WTPL WTMN ESYMWT NSYM SPINS\n", lu);
  pvinfo = vinfo;
  jj = 0;
  for (iv = 0; iv <= nvibm; ++iv)
  {
    knnmin = pvinfo->knmin[0];
    ii = pvinfo->knmin[1];
    if (ii < knnmin)
      knnmin = ii;
    ii = pvinfo->knmin[2];
    if (ii < knnmin)
      knnmin = ii;
    ii = pvinfo->knmin[3];
    if (ii < knnmin)
      knnmin = ii;
    if (knnmin < 0)
      knnmin = 0;
    k = pvinfo->gsym;
    isym = k >> 1;
    ii = (int)(pvinfo->lvupper >> 2);
    if ((pvinfo->lvupper & 2) != 0)
      ii += 10;
    if (ODD(k))
    {
      isym = -isym;
      jj = 1;
    }
    ii = ii * glob.esymdec + pvinfo->ewt[0];
    fprintf(lu, " %4d %4d %4d %4d %4d %6d %4d", iv, knnmin, pvinfo->knmax,
            pvinfo->wt[0], pvinfo->wt[3], ii, isym);
    jjs = pvinfo->spt;
    ni = jjs[0];
    for (i = 1; i < ni; ++i)
    {
      ai = jjs[i] * 0.5;
      fprintf(lu, "%5.1f", ai);
    }
    fputc('\n', lu);
    if (EVEN(isym))
    { /* even */
      i = ii;
      ii -= pvinfo->ewt[0];
      if (isym != 4)
        ii += pvinfo->ewt[1];
      if (i != ii || pvinfo->wt[0] != pvinfo->wt[2] ||
          pvinfo->wt[3] != pvinfo->wt[1])
      {
        fprintf(lu, " %4d %4d %4d %4d %4d %6d    B\n", iv, knnmin, pvinfo->knmax,
                pvinfo->wt[2], pvinfo->wt[1], ii);
      }
    }
    ++pvinfo;
  }
  if (jj != 0)
  {
    fprintf(lu, "I(TOT) IS LAST QUANTUM BEFORE F AND %s\n",
            "IS INDICATED BY NEGATIVE VALUE OF NSYM ABOVE");
  }
  /* bubble sort block number */
  lx = ntstat;
  while (lx > 1)
  {
    lv = 0;
    im = 0;
    for (i = 1; i < lx; ++i)
    {
      if (blkptr[im] > blkptr[i] ||
          (blkptr[im] == blkptr[i] && moldv[im] > moldv[i]))
      {
        lv = i;
        jj = blkptr[i];
        blkptr[i] = blkptr[im];
        blkptr[im] = jj;
        jj = moldv[i];
        moldv[i] = moldv[im];
        moldv[im] = jj;
      }
      im = i;
    }
    lx = lv;
  }
  fputs("BLOCK - WT - SYM - V - TSP - N  -  other quanta  (rel. to F=0 )\n",
        lu);
  /* convert block label to pointer */
  i = blkptr[0];
  blkptr[0] = 0;
  n = nsize = neven = nodd = maxblk = maxsblk = 0;
  for (jj = 0; jj <= ntstat; ++jj)
  {
    if (blkptr[jj] != i)
    {
      ii = jj - blkptr[n++];
      if (maxblk < ii)
        maxblk = ii;
      blkptr[n] = jj;
      if (neven > nsize)
        nsize = neven;
      if (nodd > nsize)
        nsize = nodd;
      if (jj == ntstat)
        break;
      i = blkptr[jj];
      neven = 0;
      nodd = 0;
    }
    jjt = moldv[jj];
    jjt = getqs(jjt, 0, -1, 0, ixcom, iscom, &ii);
    isym = ixcom[XSYM];
    iv = ixcom[XVIB];
    pvinfo = &vinfo[iv];
    iwt = getwt(pvinfo, isym, jjt, ivwt);
    kk = pvinfo->knmin[isym];
    if (kk < 0)
    {
      ii = ixcom[XISYM]; /* get spin symmetry */
      kk = isym;
      if (is_esym[ii] < 0)
        kk += 2;
      kk &= 2;
    }
    ii = pvinfo->knmax - kk + 2;
    if (ii < 0)
      ii = 0;
    ii >>= 1; /* size of even block */
    neven += ii;
    if (ii > maxsblk)
      maxsblk = ii;
    kk = pvinfo->knmin[3 - isym];
    if (kk < 0)
    {
      ii = ixcom[XISYM]; /* get spin symmetry */
      kk = isym;
      if (is_esym[ii] >= 0)
        kk += 2;
      kk &= 2;
    }
    ii = pvinfo->knmax - kk + 2;
    if (ii < 0)
      ii = 0;
    ii >>= 1; /* size of odd block */
    nodd += ii;
    if (ii > maxsblk)
      maxsblk = ii;
    --jjt;
    fprintf(lu, "% 5d %4d    %c %4d %4d", n + 1, iwt, csym[isym], iv, jjt);
    ni = nspin;
    for (ii = 0; ii < nspin; ++ii)
    {
      ai = iscom[ni] * 0.5;
      fprintf(lu, " %5.1f", ai);
      if (ni == itsym && itsym < itptr)
        ii = itptr;
      ni = ii;
    }
    fputc('\n', lu);
  }
  ++maxblk;
  glob.maxblk = maxblk;
  fprintf(lu, "Maximum Dimension for Hamiltonian = %d\n", nsize);
  fflush(lu);
  *nbkpf = n;
  glob.nbkpj = n;
  if (nsize < (*negy))
    *negy = nsize;
  if (ndmx > 0)
  {
    free(ivs);
    ivs = NULL;
    free(ikmin);
    ikmin = NULL;
    free(ibkptr);
    ibkptr = NULL;
    free(iqnsep);
    iqnsep = NULL;
    free(idx);
    idx = NULL;
    free(jdx);
    jdx = NULL;
    free(wk);
    wk = NULL;
    ndmx = 0;
  }
  ii = maxsblk + maxsblk - 1;
  ndmx = nsize;
  if (ndmx < ii)
    ndmx = ii;
  nl = (size_t)(nsize + ndmx) * sizeof(double);
  wk = (double *)mallocq(nl);
  wk[0] = zero;
  /* allocate space for sub-block information */
  nl = (size_t)(ndmx + 1) * sizeof(short);
  idx = (short *)mallocq(nl);
  idx[0] = 0;
  jdx = (short *)mallocq(nl);
  jdx[0] = 0;
  nl = nsize * sizeof(short);
  iqnsep = (short *)mallocq(nl);
  iqnsep[0] = 0;
  maxblk = maxblk + maxblk;
  nl = maxblk * sizeof(short);
  ibkptr = (short *)mallocq(nl);
  ibkptr[0] = 0;
  ikmin = (short *)mallocq(nl);
  ikmin[0] = 0;
  nl = (size_t)maxblk * sizeof(int);
  ivs = (int *)mallocq(nl);
  ivs[0] = 0;
  ii = glob.vibfac + 1;
  return (ii * ii);
} /* setblk */

/**
 * @brief Initializes spin factors and sorts/pre-processes Hamiltonian parameters.
 *        This function parses BCD parameter identifiers, determines operator types,
 *        calculates spin-independent and K-independent factors (zfac), and stores
 *        processed parameter information in linked lists (spar_head) for efficient
 *        Hamiltonian construction. It also handles phase conventions for operators.
 * @param list_unit File pointer for listing warnings and information. Renamed lu.
 * @param num_total_params Total number of parameters. Renamed npar.
 * @param bcd_parameter_ids Array of BCD parameter identifiers. Renamed idpar.
 * @param parameter_values Array of parameter values. Renamed par.
 * @return int Returns 1 if any diagonal K-changing operator was found (iret=1), 0 otherwise.
 *         A return of 0 might imply that glob.idiag should be set to -1 (no diagonalization).
 */
int pasort(lu, npar, idpar, par)
FILE *lu;
const int npar;
bcd_t *idpar;
const double *par;
{
  /*     initialize spin factors and arrange operators */
  /*     on entry: */
  /*         NSIZE= block size and dimension */
  /*         IDPAR= list of parameter identifiers ( element 0 is length ) */
  static char *strej[] = {"\n", "(1) vib symmetry does not match\n",
                          "(2) bad spins\n", "(3) bad spin 1 for SzNz\n",
                          "(4) specified K is beyond Kmax\n", "(5) bad weights for Itot\n",
                          "(6) specify prameters only for l >= 0\n",
                          "(7) matrix element MOD(Delta K - Delta l, nsym) is not 0\n",
                          "(8) problem with Itot symmetry\n",
                          "(9) matrix element couples states with differernt weights\n"};

  static int idmy[NDXCOM], initl;
  static bcd_t idunit[NDECPAR];
  SVIB *pvib1, *pvib2;
  SPAR *spar_now, *spar_last, *spar_match, *spar_free;
  short *iiv1, *iiv2;
  double dtmp;
  size_t nl;
  int ifac, iv1d, ipar, iv12q, ikqp, insp, njqt, ityi, isym, neulerp, i, lt;
  int kk, ltp, si1, si2, iv1, iv2, lv1, lv2, kdp, ldp, iv12, ikq, alpha, gsym;
  int ins, neuler, sznz, sznzp, si1p, si2p, nt, kd, ld, ivdif, ibcd, ldel;
  int iret, idflags, ifc, ifcp, ilim, kl, klp, nitot, ndecv, nbcd, imag;
  int alphap, ldelp, kavg, kavgp, notused, nimag[8];
  unsigned int ivsym;
  bcd_t *idval;
  BOOL first;

  if (glob.parinit > 0)
  {
    for (i = 0; i < glob.parinit; ++i)
    {
      spar_free = spar_head[i];
      while (spar_free != NULL)
      {
        spar_now = spar_free->next;
        free(spar_free);
        spar_free = spar_now;
      }
      spar_head[i] = spar_free;
    }
    spar_free = spar_now = NULL;
    glob.parinit = 0;
  }
  if (npar <= 0)
  {
    free(ipder);
    ipder = &zipder;
    initl = 0;
    return 0;
  }
  /*  find reduced spin matrix elements */
  spfac[0] = 0.;
  spfac2[0] = 0.;
  spfac[1] = sqrt(1.5);
  spfac2[0] = 0.;
  nt = glob.mxspin;
  for (i = 2; i <= nt; ++i)
  {
    dtmp = 0.5 * i; /* dtmp = I */
    spfac[i] = sqrt(dtmp * (dtmp + 1.) * (i + 1));
    spfac2[i] = 0.25 * sqrt((dtmp + 1.5) / (dtmp - 0.5)) / dtmp;
  }
  /* initialize SPECFC and DIRCOS */
  specfc(0, 0, 0, 0, 0, 0, 0, &zero, &szero, &szero);
  kk = 0;
  dircos(idmy, idmy, 0, 0, 0, &zero, &szero, &szero, 0, MODD, 0, &kk);
  ifac = glob.vibfac + 1;
  ndecv = glob.vibdec;
  ilim = ifac * ifac - 1;
  if (initl > 0)
  {
    free(ipder);
    ipder = NULL;
  }
  nl = (size_t)npar * sizeof(int);
  ipder = (int *)mallocq(nl);
  initl = npar;
  kk = 1;
  ityi = 1;
  ipder[0] = 0;
  ibcd = 0;
  nbcd = (int)idpar[0] & 0x7f;
  for (i = 1; i < npar; ++i)
  {
    ibcd += nbcd;
    if (NEGBCD(idpar[ibcd]) == 0)
    {
      ipder[i] = kk++;
      ityi = i + 1;
    }
    else
    {
      ipder[i] = -ityi;
    }
  }
  for (i = 0; i < 8; ++i)
    nimag[i] = 0;
  glob.nfit = kk;
  /* .. code parameter id  for power of N*(N+1) and equivalence */
  /* .. NJQ = POWER + 1  (first occurrance of power) */
  /* .. IPTYP= ITY + 1    (first occurrance of cosine) */
  /* ..      ITY = sequence number for hybrid operators */
  /* .. negate NJQ and IPTYP for successive occurrance */
  /* .. IP points to parameter */
  glob.parinit = glob.nvib;
  for (i = 0; i < glob.nvib; ++i)
  {
    spar_head[i] = NULL;
  }
  nsqmax = ityi = iv1d = iret = iv1 = iv2 = nitot = gsym = 0;
  ipar = npar - 1;
  ibcd = ipar * nbcd;
  spar_free = NULL;
  pvib1 = pvib2 = vinfo;
  while (ipar >= -1)
  {
    notused = -1;
    ityi = 0;
    first = TRUE;
    do
    { /* repeat until parameter subtypes exhausted */
      ivdif = 0;
      if (ipar < 0)
      { /* make sure there is a unit operator */
        idval = idunit;
        iv12q = ilim;
      }
      else
      {
        idval = &idpar[ibcd + 1];
        iv12q = bcd2i(idval[0]);
        if (ndecv > 1)
          iv12q += bcd2i(idval[1]) * 100;
        if (ndecv > 2)
          iv12q += bcd2i(idval[2]) * 10000;
        idval += ndecv;
        iv2 = iv12q / ifac;
        iv1 = iv12q - iv2 * ifac;
        ivdif = iv1 - iv2;
        if (ivdif < 0)
        {
          i = iv1;
          iv1 = iv2;
          iv2 = i;
          iv12q = iv2 * ifac + iv1;
        }
      }
      iv12 = iv12q;
      if (iv12q == ilim)
      {
        iv1 = iv2 = iv1d;
        ivdif = 0;
        iv12 = iv1d * (ifac + 1);
      }
      else if (iv1 >= glob.nvib)
      {
        /* dummy with v1 = 99 and v2 < 99 */
        if (iv1 == glob.vibfac)
          notused = 0;
        break;
      }
      if (idval[1] == (bcd_t)0x91 && idval[0] == (bcd_t)0)
      {
        /* special for rho */
        if (iv1 == iv2)
        {
          notused = 0;
          dtmp = par[ipar];
          specfc(0, iv1, iv1, 0, 0, 0, 1, &dtmp, &szero, &szero);
        }
        break;
      }
      else if (idval[1] == (bcd_t)0x99 && idval[0] == (bcd_t)0)
      {
        /* dummy parameter = 9900vv' */
        notused = 0;
        break;
      }
      pvib1 = &vinfo[iv1];
      pvib2 = &vinfo[iv2];
      gsym = pvib1->gsym;
      if ((short)gsym != pvib2->gsym)
      {
        notused = 1;
        break;
      }
      gsym = setgsym(gsym);
      nitot = glob.nitot;
      if (spar_free == NULL)
      { /* allocate more structures */
        spar_free = (SPAR *)mallocq(sizeof(SPAR));
      }
      spar_now = spar_free;
      spar_now->flags = 0;
      spar_now->alpha = C0;
      spar_now->mldel = C0;
      ityi = idpari(idval, ityi, spar_now);
      if (ityi == 0)
        break;
      njqt = (int)spar_now->njq;
      if (njqt > nsqmax)
        nsqmax = njqt;
      isym = (int)spar_now->ipsym;
      idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins, &si1, &si2,
             &sznz, &ifc, &alpha, &ldel, &kavg);
      iiv1 = pvib1->spt;
      iiv2 = pvib2->spt;
      if (si1 > 0 &&
          checksp(first, si1, si2, iiv1, iiv2, &spar_now->zfac) != 0)
      {
        notused = 2;
        break;
      }
      if (sznz > 0)
      { /*  setup for SzNz */
        if (checksp(first, 1, 0, iiv1, iiv2, &spar_now->zfac) != 0)
        {
          notused = 3;
          break;
        }
        spar_now->zfac *= 0.5;
      }
      first = FALSE;
      /* factor of 2 for B-C, etc */
      if (kd != 0 && isym == 0)
        spar_now->zfac *= 2.;
      /* setup flags */
      idflags = 1;
      imag = (kd > ld) ? kd : ld;
      if (ifc < 0)
      {
        ++imag;
        isym = 3 - isym;
        idflags |= MFCM;
      }
      if (kavg > pvib1->knmax && kavg > pvib2->knmax)
      {
        notused = 4;
        break;
      }
      if (si2 < 0)
        ++imag; /* commutator is odd */
      if (pvib1->wt[4] != pvib2->wt[4])
      {
        notused = 5;
        break;
      }
      lv1 = pvib1->lvqn;
      lv2 = pvib2->lvqn;
      if (iv12q == ilim && lv1 < 0)
        break;
      if (kavg > 0 && lv1 < 0)
      {
        if (lv1 != lv2)
          break;
        ++imag;
        isym = 3 - isym; /* 8x type operator */
      }
      if (glob.nofc && ifc != 0)
      {
        if (ifc < 0)
          ifc = -1 - ifc;
        kavg = ifc;
        ifc = 0;
        spar_now->fc = C0;
        spar_now->kavg = (char)kavg;
      }
      if (ODD(pvib1->lvupper))
      {
        --iv12;
        --iv1;
      }
      if (ODD(pvib2->lvupper))
      {
        iv12 -= ifac;
        --iv2;
      }
      /* setup for l-doubling */
      ldel = lv1 - lv2;
      if (ldel != 0)
      {
        if (ldel < 0)
          ldel = -ldel;
        if (ivdif < 0) /* DELTA K * DELTA L < 0 */
          ldel = -ldel;
        glob.lsym = FALSE;
      }
      else if (lv1 != 0 && (idflags & MLZ) == 0)
      {
        if (glob.newlz)
        {
          if (lv1 < 0)
          {
            /*  change symmetry for Lz operator */
            isym = 3 - isym;
            ++imag;
            idflags |= MLZ;
          }
        }
        else
        {
          if (lv1 < 0)
          {
            notused = 6;
            break;
          }
          i = (pvib1->lvupper ^ pvib2->lvupper) & 2;
          if (ODD2(isym + i))
          {
            /*  change symmetry for Lz operator */
            isym = 3 - isym;
            ++imag;
            idflags |= MLZ;
          }
        }
      }
      spar_now->mldel = (signed char)ldel;
      if (ODD(imag))
        idflags |= MODD;
      /* check for imaginary operators */
      if (iv1 == iv2 && ldel == 0 && isym != 0)
      {
        if (ODD(imag))
        {
          nimag[isym] += 1;
          idflags |= MIMAG;
        }
        else
        {
          nimag[isym + 4] += 1;
          idflags |= MIMAG;
        }
      }
      if (nitot < 3)
      { /* calculate alpha */
        if (gsym >= 3 && MOD(kd - ldel, gsym) != 0)
        {
          if (notused != 0)
            notused = 7;
          continue; /* matrix element breaks symmetry */
        }
        alpha = 0;
      }
      else
      {
        i = notused;
        if (notused != 0)
          notused = 8;
        alpha = MOD(ldel - kd, nitot);
        if (si1 <= itsym)
        {
          if (alpha != 0)
            continue;
          idflags |= MIDEN; /* identity spin operator */
        }
        if (alpha < 0)
          alpha += nitot;
        if (alpha == 0 || (alpha + alpha) == nitot)
        {
          if (ityi > 160)
            continue;
          if (ityi > 80 && kd == 0 && lv1 == 0 && lv2 == 0)
            continue;
        }
        else
        {
          glob.esym = FALSE;
          spar_now->zfac *= 0.5;
        }
        notused = i;
        if (ityi >= 80)
        {
          i = ityi / 80; /* i = 3,2,1,0 */
          if (ODD2(i))
          {
            if (ODD(isym + imag))
              spar_now->zfac = -spar_now->zfac;
            idflags |= MSYM2;
            isym = 3 - isym;
          }
          if (ODD(i))
          {
            alpha += nitot; /* special for quasi-diagonal */
          }
        }
      }
      if (ikq != 0 || neuler != 0 || sznz != 0 || ifc != 0 || kavg > 0)
        idflags |= MNOUNIT; /* no unit matrix */
      spar_now->alpha = (unsigned char)alpha;
      if (isym == 0 && iv1 == iv2 && kd != 0)
        iret = 1;
      if (testwt(pvib1, pvib2, isym, alpha))
      {
        if (notused != 0)
          notused = 9;
        continue;
      }
      if (ODD(neuler))
        isym = 3;
      ivsym = (unsigned int)isym + ((unsigned int)iv12 << 2);
      spar_now->ipsym = ivsym;
      spar_match = spar_last = NULL;
      for (spar_now = spar_head[iv2]; spar_now != NULL;
           spar_now = spar_now->next)
      {
        if (spar_now->ipsym <= ivsym)
          break;
        spar_last = spar_now;
      }
      if (EVEN(neuler))
      { /* not Euler denominator parameter */
        kl = idflags & MMASK;
        for (/* no init */; spar_now != NULL; spar_now = spar_now->next)
        {
          if (spar_now->ipsym < ivsym)
            break;
          spar_last = spar_now;
          klp = idpars(spar_now, &ikqp, &neulerp, &ltp, &ldp, &kdp, &insp,
                       &si1p, &si2p, &sznzp, &ifcp, &alphap, &ldelp, &kavgp);
          klp &= MMASK;
          /* check for same K dependence */
          if (ikq == ikqp && ld == ldp && kd == kdp && kl == klp &&
              sznz == sznzp && ifc == ifcp && alpha == alphap &&
              ldel == ldelp && neuler == 0 && neulerp == 0 && kavg == kavgp)
          {
            /*  check for same operator except power of N*(N+1) */
            if (lt == ltp && ins == insp && si1 == si1p && si2 == si2p)
            {
              idflags |= MNSQ;
            }
            else
            {
              if (TEST(idflags, MNSQ))
                break;
            }
            spar_match = spar_now;
          }
          else
          { /* operator with different K dependence */
            if (TEST(idflags, MNSQ))
              break;
          }
        } /* end loop over previous parameters */
      }
      notused = 0;
      if (spar_match != NULL)
      {
        if (ipar < 0)
          break; /* already have a unit operator */
        idflags |= MCOS_OK;
        spar_last = spar_match;
      }
      spar_now = spar_free;
      spar_free = NULL;
      if (spar_last != NULL)
      {
        spar_now->next = spar_last->next;
        spar_last->next = spar_now;
      }
      else
      {
        spar_now->next = spar_head[iv2];
        spar_head[iv2] = spar_now;
      }
      spar_last = NULL;
      spar_now->ip = ipar;
      spar_now->flags = (short)idflags;
    } while (ityi > 1);
    i = ipar;
    if (iv12q == ilim)
    {
      ++iv1d;
      if (iv1d >= glob.nvib)
        --ipar;
    }
    else
    {
      --ipar;
    }
    if (i > ipar)
    {
      if (notused > 0)
      {
        putbcd(sbcd, NSBCD, &idpar[ibcd]);
        if (notused * sizeof(char *) >= sizeof(strej))
          notused = 0;
        fprintf(lu,
                " WARNING: parameter %6d %s has no matrix elements: %s",
                (i + 1), sbcd, strej[notused]);
      }
      iv1d = 0;
      ibcd -= nbcd;
      notused = 0;
      ;
    }
  } /* end ipar loop */
  if (spar_free != NULL)
    free(spar_free);
  if (glob.esym)
    glob.lsym = FALSE;
  if (glob.stdphase != 0)
  { /* preset phase */
    glob.phasemask = 7;
    glob.stdphase &= 7;
  }
  else
  {
    glob.phasemask = 0;
    for (isym = 1; isym <= 3; ++isym)
    {
      /* check for phase change */
      if (nimag[isym] == 0 && nimag[isym + 4] == 0)
        continue;
      i = ipwr2[isym];
      glob.phasemask |= i;
      if (nimag[isym] > nimag[isym + 4])
      {
        glob.stdphase ^= i;
        i = nimag[isym];
        nimag[isym] = nimag[isym + 4];
        nimag[isym + 4] = i;
      }
    }
  }
  kk = glob.stdphase;
  if (kk > 0)
  {
    /* update phases */
    kk = kk ^ 5;
    for (isym = 1; isym <= 3; ++isym)
    {
      ixphase[isym] = kk & 1;
      kk = kk >> 1;
    }
  }
  nt = nimag[1] + nimag[2] + nimag[3];
  if (nt == 0)
    return iret;
  for (iv2 = 0; iv2 < glob.nvib; ++iv2)
  {
    spar_now = NULL;
    for (;;)
    {
      spar_last = spar_now;
      spar_now = (spar_now == NULL) ? spar_head[iv2] : spar_now->next;
      if (spar_now == NULL)
        break;
      isym = (int)(spar_now->ipsym & 7);
      if (isym == 0)
        continue;
      idflags = spar_now->flags;
      if ((idflags & MIMAG) == 0)
        continue;
      kk = 0;
      if (TEST(idflags, MODD))
        ++kk;
      if (TEST(glob.stdphase, ipwr2[isym]))
        ++kk;
      if (ODD(kk))
      { /* bad phase */
        ipar = spar_now->ip;
        ibcd = ipar * nbcd;
        putbcd(sbcd, NSBCD, &idpar[ibcd]);
        fprintf(lu,
                " WARNING: parameter %6d %s is imaginary and will not be used\n",
                (ipar + 1), sbcd);
        spar_free = spar_now->next;
        if (spar_last == NULL)
        {
          free(spar_now);
          spar_head[iv2] = spar_free;
        }
        else
        {
          free(spar_last->next);
          spar_last->next = spar_free;
        }
        spar_free = NULL;
        spar_now = spar_last;
        if (--nt <= 0)
          return iret;
      }
    }
  }
  return iret;
} /* pasort */

/**
 * @brief Parses the non-vibrational part of a BCD parameter identifier into a SPAR structure.
 *        It decodes fields for K^2 power, operator type, spin involvement, Fourier coefficients, etc.
 *        This function can be called multiple times for a single BCD ID if it represents
 *        multiple sub-operator types (e.g. for hybrid operators like P_a P_b + P_b P_a).
 * @param idval_bcd_param_id_no_vib Pointer to the BCD array for the non-vibrational part of the parameter ID. This array is NDECPAR long.
 * @param itp_subtype_in Input sub-parameter type. If 0, it's the first call for this BCD ID.
 *                       If >0, it's a subsequent call to get the next sub-operator component.
 * @param pspar_out Output: Pointer to the SPAR structure to be filled with parsed fields.
 * @return int The next itp_subtype_in value to use for further calls (if more sub-operators exist for this ID),
 *             or 0 if no more sub-operators for this BCD ID or if an invalid combination is encountered.
 */
int idpari(idval, itp, pspar)
bcd_t *idval;
int itp;
SPAR *pspar;
{
  /* subroutine to parse parameter identifier for complex parameter type */
  /*     on entry: */
  /*         IDVAL= parameter ID / vibrational field */
  /*         ITP= value of sub-parameter type for hybrid operator */
  /*     on return: */
  /*         pspar modified */

  /*0  1  2   3   4   5   6   7  8   9 */
  static int itpnxt[10] = {0, 0, 0, 2, 6, 6, 2, 0, 9, 0};
  static int itpop[10] = {0, 0, 0, 10, 40, 40, 10, 1, 0, 1};
  static int zfacv[10] = {1, 2, 6, 3, -8, 8, -6, 2, 4, -4};
  int ins, ity, nsx, isy, itysav;
  int itmp, i, iphaz, ldt, lnt, kdt, ksq, si1, si2, exval;
  bcd_t ibtmp;

  iphaz = 0;
  if (itp > 0)
  {
    /* get current sub-parameter type from previous */
    if (itp < 10)
    {
      itp = itpnxt[itp];
      if (itp == 0)
        return 0;
    }
    else
    {
      iphaz = itp >> 4;
      itp = itpnxt[itp & 15];
    }
  }
  pspar->zfac = 1.;
  pspar->euler = C0;
  pspar->msznz = C0;
  pspar->fc = C0;
  exval = 0;
  ibtmp = idval[0];
  nsx = (int)(ibtmp & 0x0f);
  ksq = (int)(ibtmp >> 4);
  ity = bcd2i(idval[1]);
  itysav = ity * 10 + ksq;
  ibtmp = idval[2];
  ins = (int)(ibtmp & 0x0f);
  si1 = (int)(ibtmp >> 4);
  if (ins > 0 && nspin == 0)
    return 0;
  if (ins >= 5 && itsym == 0)
    return 0;
  ibtmp = idval[3];
  si2 = (int)(ibtmp & 0x0f);
  itmp = (int)(ibtmp >> 4);
  if (si2 > si1)
  {
    /*  make sure SI1 > SI2 */
    i = si1;
    si1 = si2;
    si2 = i;
    if (si2 == 0)
      si2 = -1;
  }
  if (si1 > nspin)
    return 0;
  if (ity > 90)
  {
    if (si1 != 0 || ins >= 5)
      return 0;
  }
  ibtmp = idval[4];
  if (ibtmp != (bcd_t)0)
  {
    if (ibtmp >= (bcd_t)0x60)
      return 0;
    itmp += 10 * (int)(ibtmp & 0x0f);
    exval = (int)ibtmp & 0xf0;
  }
  pspar->kavg = C0;
  if (itmp != 0)
  {
    i = (itmp - 1) / 10;
    itmp -= (i >> 1) * 10;
    if (ODD(i))
      itmp = 9 - itmp;
    pspar->fc = (char)itmp;
  }
  pspar->msi1 = (char)si1;
  pspar->msi2 = (char)si2;
  pspar->njq = (unsigned char)nsx;
  /*        exclude S * S and I * S from aa, bb, cc operators */
  if (itp == 2)
  {
    if (si2 > 0)
    {
      itp = itpnxt[itp];
      if (itp == 0 && iphaz == 0)
        return 0;
    }
  }
  if (itp != 0)
  {
    if (itp == 6 && si1 == 0)
    {
      itp = 8;
      pspar->njq = (unsigned char)(nsx + 1);
    }
  }
  else
  {
    if (iphaz != 0)
    {
      if (ins >= 5)
        iphaz -= 1;
      else if (glob.nitot >= 3)
        iphaz -= 5;
    }
    else
    {
      if (ins >= 5)
        iphaz = 4; /* iphaz = 4,3,2,1,0 */
      if (glob.nitot >= 3)
        iphaz += 15;
    }
    if (ity == 0)
    {
      itp = 1;
    }
    else if (ity <= 3)
    {
      itp = ity;
      if (glob.oblate)
        itp = revsym[itp];
      itp += 2;
      /*  ITP=3,4,5 for A,B,C */
      if (itp == 3 && si1 == 0)
        itp = 7;
    }
    else
    {
      itp = 1;
      itmp = ity;
      if (ity >= 92)
        itmp = (ity - 92) >> 1;
      /*  negate odd powers of P+**2 */
      if (itmp <= 11 && EVEN(itmp))
        pspar->zfac = -(pspar->zfac);
    }
  }
  isy = ity / 20;
  if (isy > 3)
  {
    isy = 0;
  }
  else if (isy > 0)
  {
    if (glob.oblate)
    {
      itysav += (revsym[isy] - isy) * 200;
    }
    else
    {
      isy = revsym[isy];
    }
  }
  pspar->ipsym = (unsigned int)isy;
  if (itp > 1)
  {
    itysav += itpop[itp] - ity * 10;
    pspar->zfac *= 2. / (double)zfacv[itp];
  }
  itp += iphaz << 4;
  if (ins >= 5)
  {
    ins -= 5;
    if (iphaz >= 5)
      iphaz = MOD(iphaz, 5);
    pspar->msznz = (unsigned char)(iphaz + 1);
  }
  pspar->mins = (unsigned char)ins;
  ity = itysav / 10;
  ksq = itysav - ity * 10;
  ldt = kdt = 0;
  if (ity <= 3)
  {
    if (ity != 0)
      ldt = 2;
  }
  else if (ity < 12)
  {
    ldt = kdt = ity + ity - 6;
  }
  else if (ity < 20)
  {
    kdt = ity + ity - 22;
    ldt = kdt + 1;
  }
  else if (ity < 80)
  {
    isy = ity / 20;
    ldt = kdt = ity - isy * 20 + 1;
    if (ODD((isy >> 1) + ldt))
      --kdt;
  }
  else if (ity < 90)
  {
    /* K energies coded as 8n  with K=10*n+KSQ */
    if (exval != 0)
      return 0;
    pspar->kavg = (char)(10 * (ity - 80) + ksq);
    ksq = 0;
  }
  else
  {
    if (exval != 0)
      return 0;
    i = ity - 90;
    if (ODD(i))
    {
      if (ksq > 1)
      { /* augmented Euler series */
        --i;
        ksq += 8;
        kdt = i;
      }
    }
    else
    {
      kdt = i;
    }
    pspar->euler = (unsigned char)(i + 2);
  }
  pspar->ksq = (char)ksq;
  if (si1 == 0)
  {
    lnt = 0;
  }
  else if (si2 <= 0)
  {
    if (isy == 0 || ldt == 0)
      si2 = 0;
    else if (si2 < 0)
      --ldt;
    lnt = 1;
  }
  else if (ldt < 2)
  {
    lnt = 0;
  }
  else
  {
    lnt = 2;
  }
  pspar->mln = (unsigned char)lnt;
  pspar->mkdel = (unsigned char)kdt;
  /*  try to use raising operator */
  if (ldt == kdt && ldt > lnt)
    ldt -= (ldt - lnt) & 0x7e; /* subtract even number */
  pspar->mld = (unsigned char)ldt;
  if (exval != 0)
  {
    /* check for Euler restrictions */
    ibtmp = (bcd_t)exval;
    exval = exval >> 3;
    if (idval[1] == (bcd_t)0 && idval[2] == (bcd_t)0 &&
        idval[3] == (bcd_t)0 && idval[4] == ibtmp)
    {
      /* energy Euler term */
      if (idval[0] == (bcd_t)0)
        return 0; /* bad energy-like term */
    }
    else if (exval == 2)
    {
      exval = 12; /* operator Euler term */
    }
    pspar->euler = (unsigned char)exval;
  }
  return itp;
} /* idpari */

/**
 * @brief Checks spin selection rules and calculates part of the spin-dependent factor.
 *        Also initializes I_tot related coefficients via setzitot for the first call for a spin set.
 * @param first_call_for_spin_set TRUE if this is the first time checksp is called for this specific
 *                                combination of spin indices (si1, si2) and vibrational states.
 *                                Used to initialize I_tot factors. Renamed first.
 * @param spin1_idx Index of the first spin involved in the operator.
 * @param spin2_idx Index of the second spin involved (0 if one-spin, >0 if two-spin).
 * @param vib_bra_spin_table Spin table (spt) for the bra vibrational state.
 * @param vib_ket_spin_table Spin table (spt) for the ket vibrational state.
 * @param matrix_element_factor_ptr Pointer to the matrix element factor, which will be multiplied by spin-related prefactors. Renamed zfac.
 * @return int 0 if spin selection rules are satisfied, otherwise a non-zero error code (2-6, or 1 for I_tot related issue).
 */
int checksp(const BOOL first, int si1, int si2, const short *iiv1,
            const short *iiv2, double *zfac)
{
  int ii, iip;
  ii = iip = 0;
  if (si2 > 0)
  {
    if ((short)si2 < iiv1[0])
      ii = iiv1[si2];
    if ((short)si2 < iiv2[0])
      iip = iiv2[si2];
    if (ii == iip)
    {
      if (ii == 0)
        return 2;
      if (si2 == si1)
      {
        if (ii < 2)
          return 3;
        *zfac *= spfac2[ii];
      }
      *zfac *= spfac[ii];
    }
    else
    { /* ii != iip */
      if (ODD(ii + iip))
        return 4; /* check multiplicity */
      if (si2 > itsym)
        return 4;
    }
  }
  if (si2 != si1)
  {
    ii = iip = 0;
    if ((short)si1 < iiv1[0])
      ii = iiv1[si1];
    if ((short)si1 < iiv2[0])
      iip = iiv2[si1];
    if (ii == iip)
    {
      if (ii == 0)
        return 5;
      *zfac *= spfac[ii];
    }
    else
    { /* ii != iip */
      if (ODD(ii + iip))
        return 6; /* check multiplicity */
      if (si1 > itsym)
        return 6;
    }
  }
  if (first && si1 > itsym)
  {
    iip = 0;
    if (glob.nitot >= 3)
      iip = si1 - itsym;
    if (si1 == si2)
    {
      if (iip > 1)
        return 1;
      setzitot(2, 0, 2, ii, glob.nitot); /* quadrupole */
    }
    else if (si2 > itsym)
    {
      if (iip > 2)
        return 1;
      setzitot(1, 1, 0, ii, glob.nitot); /* 2-spin product */
      setzitot(1, 1, 2, ii, glob.nitot);
    }
    else
    {
      if (iip > 1)
        return 1;
      setzitot(1, 0, 1, ii, glob.nitot); /* 1-spin vector */
    }
  }
  return 0;
} /* checksp */

/**
 * @brief Sets statistical weights for a vibrational state based on symmetry options.
 * @param pvinfo_vib_state Pointer to the SVIB structure for the vibrational state.
 * @param vib_state_index_or_count If positive, it's the current vibrational state index.
 *                                 If negative, it's -N_total_vib_states, indicating to apply
 *                                 these weights to N_total_vib_states starting from pvinfo_vib_state.
 * @param axis_option_iax Internal axis code (1=z(a), 2=y(b), 3=x(c), 4=A_D2h_like, 5=B_D2h_like) for D2 or Cn symmetries.
 * @param weight_plus Statistical weight for "plus" or "even" parity states.
 * @param weight_minus Statistical weight for "minus" or "odd" parity states.
 * @param vibrational_symmetry_code_vsym Encodes vibrational symmetries for multiple states if vib_state_index_or_count is negative,
 *                                       or indicates if current state's weights should be swapped if positive.
 * @return int Always 0.
 */
int setwt(pvinfov, ivib, iax, iwtpl, iwtmn, vsym)
SVIB *pvinfov;
const int ivib, iax, iwtpl, iwtmn;
double vsym;
{
  /*  set weights */
  /*     on entry: */
  /*         IVIB= number of vibrational states * number of wang states */
  /*         IAX= axis of symmetry */
  /*         IWTPL= weight for even rotational states */
  /*         IWTMN= weight for odd rotational states */
  /*         VSYM= symmetry code for vibrational states */
  /*     on return: */
  /*         IWT= weights for vibration rotation states */
  SVIB *pvinfo;
  int n, nvsym, i, k;
  short ii, jj, itmp;

  n = 1;
  nvsym = 0;
  if (ivib < 0)
  {
    if (vsym > 0.5)
      nvsym = 1;
    n = -ivib;
  }
  pvinfo = pvinfov;
  for (i = 0; i < n; ++i)
  {
    ii = (short)iwtpl;
    jj = (short)iwtmn;
    if (nvsym != 0)
    {
      k = (int)(fmod(vsym, 10.) + 0.5);
      vsym = (vsym - k) / 10.;
      if (vsym < 0.)
        nvsym = 0;
      if (ODD(k))
      {
        itmp = ii;
        ii = jj;
        jj = itmp;
      }
    }
    switch (iax)
    {
    case 1: /*  WEIGHTS FOR Z */
      pvinfo->wt[0] = ii;
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = jj;
      pvinfo->wt[3] = ii;
      pvinfo->wt[4] = 2;
      ii = jj;
      break;
    case 2: /* WEIGHTS FOR B */
      pvinfo->wt[0] = ii;
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = ii;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 3;
      break;
    case 3: /* WEIGHTS FOR X */
      pvinfo->wt[0] = ii;
      pvinfo->wt[1] = ii;
      pvinfo->wt[2] = jj;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 2;
      break;
    case 4: /* A1, A2 only */
      pvinfo->wt[0] = ii;
      pvinfo->wt[3] = jj;
      pvinfo->wt[4] = 3;
      break;
    case 5: /* B1, B2 only */
      pvinfo->wt[1] = jj;
      pvinfo->wt[2] = ii;
      pvinfo->wt[4] = 3;
      break;
    } /* end switch */
    ++pvinfo;
  }
  return 0;
} /* setwt */

/**
 * @brief Gets the appropriate statistical weight for a given state symmetry and spin configuration.
 * @param pvinfo_vib_state Pointer to the SVIB structure for the vibrational state.
 * @param d2_symmetry_idx D2 symmetry index (0-3) of the rotational state.
 * @param spin_state_idx If 0, requests general check for existence of states for this d2_symmetry_idx.
 *                       If >0, it's the spin pattern index (1-based) from the spin table (spt).
 *                       If <0, it's -alpha_symmetry_code for I_tot E-state weight lookup.
 * @param weights_out Output array of 3 integers:
 *                    weights_out[0] = weight for A/B type rotational states (or primary component).
 *                    weights_out[1] = weight for E type rotational states (first E component, e.g. E_a).
 *                    weights_out[2] = weight for other E or B type states (e.g. E_b or B for D4/D6).
 *                    Unused components are set to -1.
 * @return int The primary statistical weight (weights_out[0]), or sum of weights if I_tot basis and multiple components allowed.
 *             Returns <=0 if no states/weight for this symmetry/spin.
 */
int getwt(pvinfo, isym, iispin, ivwt)
SVIB *pvinfo;
const int isym, iispin;
int *ivwt;
{
  int iwt, jsym, nset, k, msym, nsym;
  short *jjs;
  if (iispin == 0)
  {
    iwt = 0;
    if (pvinfo->knmax >= pvinfo->knmin[isym] ||
        pvinfo->knmax >= pvinfo->knmin[3 - isym])
      iwt = 1;
    return iwt;
  }
  k = 0;
  jsym = isym & 3;
  iwt = pvinfo->wt[4];
  if (iwt != 0)
  {
    /* resolve Itot symmetries */
    if (iispin > 0)
    {
      jjs = pvinfo->spt;
      nset = jjs[0];
      k = jjs[iispin * nset + itsym + 1];
    }
    else
    {
      k = -iispin;
    }
    if (iwt > 0 && ODD2(k))
      jsym ^= iwt;
  }
  ivwt[0] = pvinfo->wt[jsym];
  ivwt[1] = ivwt[2] = -1;
  nsym = pvinfo->gsym >> 1;
  if (nsym < 3)
    return ivwt[0];
  if (ODD(nsym))
  {                           /* nsym = 3, 5 */
    ivwt[1] = pvinfo->ewt[0]; /* E state */
  }
  else
  { /* nsym = 4, 6 */
    msym = isoddk[jsym] - (int)(pvinfo->lvqn);
    if (iwt != 0)
      msym += k >> 2;
    if (ODD(msym))
    {
      ivwt[0] = -1;             /* no A state */
      ivwt[1] = pvinfo->ewt[1]; /* E state */
      if (nsym == 6)
        ivwt[2] = pvinfo->wt[jsym]; /* B state */
    }
    else
    {
      ivwt[1] = pvinfo->ewt[0]; /* E state */
      if (nsym == 4)
        ivwt[2] = pvinfo->wt[jsym ^ 2]; /* B state */
    }
  }
  iwt = ivwt[0];
  if (pvinfo->gsym > 6)
  {
    if (iwt <= 0)
      iwt = ivwt[1];
    if (iwt <= 0)
      iwt = ivwt[2];
  }
  return iwt;
} /* getwt */

/**
 * @brief Tests if an operator can connect two vibrational states based on their statistical weights
 *        and the operator's symmetry.
 * @param pvib1_bra Pointer to SVIB structure for the bra state.
 * @param pvib2_ket Pointer to SVIB structure for the ket state.
 * @param operator_d2_symmetry D2 symmetry of the operator (0-3).
 * @param operator_alpha_symmetry I_tot alpha symmetry component of the operator.
 * @return BOOL TRUE if the connection is forbidden by weight/symmetry rules, FALSE otherwise.
 */
BOOL testwt(pvib1, pvib2, isym, alpha)
SVIB *pvib1, *pvib2;
int isym, alpha;
{ /* return true if all weight sets do not match */
  int ii, jj, kk, ix, jx, iwt4, iwt[3], jwt[3];
  if (pvib1->gsym != pvib2->gsym)
    return TRUE;
  iwt4 = pvib1->wt[4];
  if ((short)iwt4 != pvib2->wt[4])
    return TRUE;
  ix = -1;
  jx = ix - (alpha << 2);
  for (kk = 0; kk < 2; ++kk)
  {
    for (ii = 0; ii < 4; ++ii)
    {
      getwt(pvib1, ii, ix, iwt);
      jj = ii ^ isym;
      getwt(pvib2, jj, jx, jwt);
      if (checkwt(iwt, jwt) == 0)
        return FALSE;
    }
    if (iwt4 <= 0)
      break;
    jx -= 2;
  }
  return TRUE;
} /* testwt */

/**
 * @brief Checks if any component of two sets of statistical weights are equal and positive.
 *        Used by testwt to determine if an operator can connect states.
 * @param weights_bra Array of 3 weights for the bra state (A/B, E1, E2/B').
 * @param weights_ket Array of 3 weights for the ket state.
 * @return int 0 if any pair of corresponding positive weights are equal, 1 otherwise.
 */
int checkwt(iwt, jwt)
int *iwt, *jwt;
{
  /* return 0 if any pair of weights are equal */
  int ii, jj, nn;
  nn = glob.maxwt;
  if (nn == 0)
    return (iwt[0] - jwt[0]);
  for (ii = nn; ii >= 0; --ii)
  {
    if (iwt[ii] <= 0)
      continue;
    for (jj = nn; jj >= 0; --jj)
    {
      if (iwt[ii] == jwt[jj])
        return 0;
    }
  }
  return 1;
} /* checkwt */

/**
 * @brief Sets the global symmetry context (n-foldness, I_tot parameters) based on a gsym code.
 * @param group_symmetry_code The gsym code for a vibrational state (from SVIB struct).
 *                            Even values are standard symmetry, odd values indicate I_tot basis active.
 *                            Magnitude >> 1 gives n-foldness.
 * @return int The n-foldness of the symmetry group.
 */
int setgsym(const int gsym)
{
  static int oldgsym = -1;
  static int nsym;
  int k, kk;
  if (gsym == oldgsym)
    return nsym;
  oldgsym = gsym;
  nsym = gsym >> 1;
  is_esym[0] = 0;
  for (k = 1, kk = nsym - 1; k < kk; ++k, --kk)
  {
    is_esym[k] = 1;
    is_esym[kk] = -1;
  }
  if (k == kk)
    is_esym[k] = 0;
  if (ODD(gsym))
  {
    itptr = nspin - 2;
    itsym = nspin - nsym;
    glob.nitot = nsym;
  }
  else
  {
    itptr = itsym = nspin + nspin + 1;
    glob.nitot = 0;
  }
  /* set up nominal spin coupling map */
  for (k = 0; k < nspin; ++k)
  {
    ismap[k + k] = k - 1; /* last N,J,F1 .. */
    if (k < itsym)
    {
      ismap[k + k + 1] = k + nspin + 1; /* S,I1,I2 .. */
    }
    else
    {
      ismap[k + k + 1] = itptr; /* Itot */
      break;
    }
  }
  ismap[0] = nspin; /* fix up position of N */
  return nsym;
} /* setgsym */

/**
 * @brief Parses BCD encoded spin quantum numbers and associates them with a vibrational state.
 *        It manages a list of unique spin patterns (SSP structs) to save memory.
 * @param bcd_spin_codes Pointer to the BCD array containing spin codes. The first element [0]
 *                       holds flags and length. Subsequent elements encode 2*I values. Renamed ispnx.
 * @param pvinfo_vib_state Pointer to the SVIB structure for the vibrational state to associate
 *                         the spins with. Renamed pvinfo.
 * @return int Number of actual spin quantum numbers defined for this state (excluding N),
 *             adjusted for I_tot symmetry if applicable.
 */
int getsp(ispnx, pvinfo)
    const bcd_t *ispnx;
SVIB *pvinfo;
{
  /*  subroutine to set up spin structure                        */
  /*     on entry:                                               */
  /*         ISPNX = spin code                                   */
  /*         IVIB  = vibrational quantum number                  */
  /*     on return:                                              */
  /*         JJS array of spin combinations                      */
  /*         (this array is initialized mostly in setsp)         */
  /*                                                             */
  SSP *ssp_now;
  short *jjs;
  size_t nl;
  int i, ii, jj, nspinv, idec, nspnx, nsstat, icmp, nset, isbig;
  int nitot, iret;
  short iis[MAXSPIN + 1];
  bcd_t itmp;

  jj = -1;
  nspnx = (int)(ispnx[0] & 0x7f);
  isbig = 0;
  itmp = ispnx[1];
  if ((int)(itmp & 0x0f) == 0)
  {
    jj = (int)(itmp >> 4);
    isbig = 1;
  }
  /*  decode spins */
  nspinv = 0;
  nl = 1;
  idec = isbig;
  for (i = 1; i <= MAXSPIN; ++i)
  {
    if (isbig == 0 && jj >= 0)
    {
      ii = jj;
      jj = -1;
    }
    else
    {
      ++idec;
      if (idec >= nspnx)
        break;
      itmp = ispnx[idec];
      ii = (int)(itmp & 0x0f);
      if (isbig != 0)
        ii = ii * 10 + jj;
      jj = (int)(itmp >> 4);
    }
    if (ii == 0)
      break;
    if (ii >= MAXII)
    {
      ii = 1;
    }
    else if (ii > 1)
    {
      nspinv = i;
    }
    nl *= (size_t)ii;
    --ii;
    iis[i] = (short)ii;
    if (ii > glob.mxspin)
      glob.mxspin = ii;
  }
  nsstat = (int)nl;
  i = (int)(nsstat << glob.msshft);
  if (i < 0 || (size_t)(i >> glob.msshft) != nl)
  {
    puts("spin problem too big");
    exit(EXIT_FAILURE);
  }
  nitot = 0;
  iret = nspinv;
  i = pvinfo->gsym;
  if (ODD(i))
  {
    nitot = i >> 1;
    if (nitot > nspinv)
    {
      iret += nitot;
      if (iret > MAXSPIN)
        iret = MAXSPIN;
      for (i = nspinv + 1; i <= iret; ++i)
        iis[i] = 0;
      nspinv = iret;
    }
    if (nitot > 3)
    {
      iret -= (nitot - 3);
    }
  }
  pvinfo->nspstat = nsstat;
  nset = nspinv + 1;
  iis[0] = (short)nset;
  ssp_now = &ssp_head;
  do
  { /*  check previous combinations for match */
    icmp = 0;
    jjs = ssp_now->sspt;
    for (i = 0; i < nset; ++i)
    {
      icmp = iis[i] - jjs[i];
      if (icmp != 0)
        break;
    }
    if (icmp == 0)
      icmp = nitot - ssp_now->nitot;
    /*  match obtained, so return */
    if (icmp == 0)
    {
      pvinfo->spt = jjs;
      return iret;
    }
    ssp_now = ssp_now->next;
  } while (ssp_now != NULL);
  nl = (size_t)(nsstat + 1) * (size_t)nset * sizeof(short);
  jjs = (short *)mallocq(nl);
  nl = sizeof(SSP);
  ssp_now = (SSP *)mallocq(nl);
  ssp_now->next = ssp_head.next;
  ssp_head.next = ssp_now;
  ssp_now->sspt = jjs;
  ssp_now->ssize = nsstat;
  ssp_now->nitot = nitot;
  jjs[0] = (short)nset;
  for (i = 1; i < nset; ++i)
  {
    jjs[i] = iis[i];
  }
  pvinfo->spt = jjs;
  return iret;
} /* getsp */

/**
 * @brief Finalizes the spin tables (spt) for all unique spin patterns (SSP structs).
 *        For each unique spin pattern, it calculates and stores the coupling information:
 *        - Differences like 2*(N-F), 2*(J-F), 2*(F1-F), etc.
 *        - Minimum 2*F value for each spin combination.
 *        - I_tot symmetry components if applicable.
 *        This pre-computation is used by getqs to quickly retrieve quantum numbers.
 */
void setsp(void)
{
  /*  subroutine to set up spin structure                        */
  /*   let: jjs = array of spin combinations in struct ssp.sspt  */
  /*         jjs[0]  = NSET = NSPIN + 1 = size of the set        */
  /*                   {N,J,F1..F} or {N,J..Itot,F}              */
  /*         jjs[1]..jjs[NSPIN] = S, I1, I2..                    */
  /*         jjs[k*NSET]..jjs[NSPIN-1+k*NSET] = 2*(N-F),2*(J-F)..*/
  /*         jjs[NSPIN+k*NSET] = minimum value of 2*F +(1)       */
  /*                                                             */
  SSP *ssp_now;
  ITMIX *pitmix0;
  EITMIX *pitmix;
  int nset, nspinv, ispsum, nn, iv, itot, ibgn, iend, ns, ns0, ifhalf;
  int minf0, minft, jj, ival, ibase, i, itmp, knt, nlim, nitot, ii;
  int itsym0, itptr0;
  short *iis, *jjs, *jjs0;
  nspin = 0;
  ssp_now = &ssp_head; /* head has no spins */
  while ((ssp_now = ssp_now->next) != NULL)
  { /* loop over spin sets */
    iis = ssp_now->sspt;
    nspinv = iis[0] - 1;
    for (iv = nspinv; iv > nspin; --iv)
    {
      if (iis[0] > 0)
      {
        nspin = iv;
        break;
      }
    }
  }
  ssp_now = &ssp_head; /* head has no spins */
  while ((ssp_now = ssp_now->next) != NULL)
  { /* loop over spin sets */
    iis = ssp_now->sspt;
    nset = iis[0];
    nspinv = nset - 1;
    if (nspinv > nspin)
      nspinv = nspin;
    ispsum = 0;
    for (iv = 1; iv <= nspinv; ++iv)
    {
      ispsum += iis[iv];
    }
    ifhalf = (ispsum & 1);
    iend = iis[nspinv];
    ival = iend + 1;
    ns = ival;
    pitmix0 = NULL;
    pitmix = NULL;
    nitot = ssp_now->nitot;
    if (nitot == 0)
    {
      nlim = nspinv;
      ibgn = iend;
      itsym0 = itptr0 = nspin + nspin + 1;
    }
    else
    { /* Itot coupling */
      itptr0 = nspin - 2;
      itsym0 = nspin - nitot;
      nlim = itsym0 + 1;
      ii = iend;
      for (i = nlim; i < nspinv; ++i)
      {
        if ((int)iis[i] != ii)
          break;
        ns *= ival;
      }
      if (i < nspin)
      {
        puts("spins under Itot should be the same");
        exit(EXIT_FAILURE);
      }
      pitmix0 = get_itmix(ii, nitot);
      if (nitot < 3)
        pitmix0 = NULL;
      iend = nitot * ii;
      ibgn = (iend & 1);
    }
    knt = ssp_now->ssize;
    ns0 = knt / ns;
    jjs = iis + nset;
    for (nn = -ispsum; nn <= ispsum; nn += 2)
    { /* select N quanta */
      minf0 = -ifhalf;
      if (nn < minf0)
        minf0 = nn; /* nn = 2 * (N - F) */
      if (pitmix0 != NULL)
        pitmix = pitmix0->mix;
      for (itot = ibgn; itot <= iend; itot += 2)
      {
        ns = ns0 * (itot + 1);
        if (pitmix != NULL && itot > ibgn)
          pitmix += 1;
        for (iv = 0; iv < ns; ++iv)
        { /* try all possible spins */
          jj = nn;
          jjs[0] = (short)nn;
          ival = iv;
          minft = minf0;
          for (i = 1; i < nlim; ++i)
          {
            ibase = iis[i] + 1;
            itmp = ival;
            ival /= ibase;
            itmp -= ival * ibase;
            jj += itmp - iis[i]; /* jj = J + N - S - 2F */
            if (jj < minft)
              minft = jj;
            jj += itmp; /* jj = 2 * (J - F) */
            jjs[i] = (short)jj;
          }
          jj += ival - itot;
          if (jj < minft)
            minft = jj;
          jj += ival;
          if (jj == 0)
          { /* this possibility gives F=0 */
            if (itptr0 < nspinv)
              jjs[itptr0 + 1] = (short)itot;
            jj = ifhalf - minft;
            minft = jj + (jj & 1);
            jjs[nspinv] = (short)minft;
            if (pitmix != NULL)
            {
              jjs[itsym0 + 1] = pitmix->qsym[0];
              for (i = itsym0 + 2; i <= itptr0; ++i)
                jjs[i] = 0;
              nitot = pitmix->n;
              jjs0 = jjs;
              for (i = 1; i < nitot; ++i)
              {
                jjs += nset;
                --knt;
                for (jj = 0; jj <= nspinv; ++jj)
                  jjs[jj] = jjs0[jj];
                jjs[itsym0 + 1] = pitmix->qsym[i];
              }
            }
            if (--knt <= 0)
            {
              ns0 = 0;
              break;
            }
            jjs += nset;
          }
        } /* end loop over all spins */
      } /* end loop over Itot */
    } /* end loop over N */
  } /* end loop over spin types */
} /* setsp */

/**
 * @brief Gets the sub-block structure and key quantum numbers for a given F-block.
 * @param block_index Index of the F-block (0 to num_F_blocks-1).
 * @param f_qn_times_2_out Output: Pointer to store 2*F for this block.
 * @param stat_weights_block_out Output: Array to store statistical weights [A,E1,E2/B,gsym,nqn] for this block.
 * @param sub_block_pointers_out Output: Array to store start index of each sub-block (Wang*spin) within this F-block.
 * @param min_k_values_out Output: Array to store minimum K value for each sub-block.
 * @param vib_spin_sym_packed_out Output: Array to store packed (Vib,Sym,SpinPattern) identifier for each sub-block.
 * @return int Number of sub-blocks in this F-block. Returns 0 if block_index is invalid or no states.
 */
int getqq(iblk, f, iwtb, sbkptr, kmin, vs)
    const int iblk;
int *f, *iwtb, *vs;
short *sbkptr, *kmin;
{
  /* subroutine to get sub-block structure and major quantum numbers */
  /*     on entry: */
  /*         MXBLK=biggest block size */
  /*         IBLK= block number */
  /*     on return: */
  /*         GETQQ= number of sub-blocks */
  /*         F= F quantum number *2 */
  /*         IWTB= statistical weight */
  /*         SBKPTR= pointer to beginning index for sub-block */
  /*         KMIN= minimum K for each sub-block */
  /*         VS= vibration and symmetry for each sub-block */
  SVIB *pvinfo;
  short *jjs;
  int ibgn, iend, ksym, i, n, nsize, isblk, nsblk, k, kk, nn, iv;
  int mss, nset, nsym, ff;
  unsigned int mvs;

  isblk = iblk - 1;
  if (glob.nbkpj <= 0)
    return 0;
  ff = isblk / glob.nbkpj;
  isblk -= glob.nbkpj * ff; /* form remainder */
  ff = ff + ff;
  /*  get information for sub-blocks */
  nsize = nsblk = nsym = 0;
  ibgn = blkptr[isblk];
  iend = blkptr[isblk + 1];
  for (i = ibgn; i < iend; ++i)
  {
    mvs = (unsigned int)moldv[i];
    ksym = (int)mvs & 3;
    iv = (int)(mvs >> 2) & glob.msmask;
    mss = (int)(mvs >> glob.msshft);
    pvinfo = &vinfo[iv];
    if (i == ibgn)
      nsym = setgsym((int)pvinfo->gsym);
    jjs = pvinfo->spt;
    nset = jjs[0];
    jjs += nset * mss;
    if (ff < jjs[nset - 1])
      continue;
    nn = jjs[0] + ff;
    n = nn >> 1;
    if (ODD(n))
      ksym = 3 - ksym;
    k = pvinfo->knmin[ksym];
    if (k < 0)
    {
      kk = MOD(jjs[itsym + 1] >> 2, nsym); /* get spin symmetry */
      k = ksym;
      if (is_esym[kk] < 0)
        k += 2;
      k &= 2;
    }
    kk = pvinfo->knmax;
    if (kk > n)
      kk = n;
    kk -= k;
    if (kk >= 0)
    {
      sbkptr[nsblk] = (short)nsize;
      kmin[nsblk] = (short)k;
      vs[nsblk] = (int)mvs;
      nsize += 1 + (kk >> 1);
      if (nsblk == 0)
      {
        *f = ff - (nn & 1);
        kk = -4;
        if (pvinfo->wt[4] != 0)
          kk -= jjs[itsym + 1];
        getwt(pvinfo, (int)(mvs & 3), kk, iwtb);
        iwtb[3] = pvinfo->gsym;
        iwtb[4] = pvinfo->nqn;
      }
      ++nsblk;
    }
  }
  if (nsblk > 0)
    sbkptr[nsblk] = (short)nsize;
  return nsblk;
} /* getqq */

/**
 * @brief Gets the full set of quantum numbers for a specific state within a block.
 * @param f_block_idx F-block index. Renamed iblk.
 * @param state_idx_in_block Index of the state within the F-block (0 to Nblock-1). Renamed indx.
 * @param max_qn_fields Max number of quantum number fields to fill in qn_list_out. Renamed maxqn.
 * @param qn_list_out Output array to store the quantum numbers. Renamed iqn.
 * @param degeneracy_out Output: Pointer to store the degeneracy of the state. Renamed idgn.
 * @return int Returns a code related to Wang block sorting or K-parity (ncod).
 *             Positive for one K-parity, negative for the other, for asymmetric tops.
 *             1 for symmetric top / linear K=0, or if K-parity not distinguished.
 */
int getqn(iblk, indx, maxqn, iqn, idgn)
    const int iblk,
    indx, maxqn;
short *iqn;
int *idgn;
{
  /*   subroutine to get quantum numbers */
  /*     on entry: */
  /*         IBLK= block number */
  /*         INDX= index within block */
  /*     on return: */
  /*         IQN= list of quantum numbers */
  /*         IDGN= degeneracy */
  /*         NCOD= seraching aid to find state within a Wang block */
  static int last;
  int ibgn, iend, iqsp, isym, i, k, nsblk, ioff, lower, iv, iqf;
  int n, ncod, ldgn, ix, kq, ivbase, lv;
  short mgsym, ngsym;

  /*  get sub-block info from store or GETQQ */
  /*     look at last used then oldest */
  if (iblk == cgetv[0].cblk)
  {
    cgetq = cgetv;
    last = 0;
    ioff = 0;
    nsblk = cgetq->cnblk;
  }
  else if (iblk == cgetv[1].cblk)
  {
    cgetq = &cgetv[1];
    last = 1;
    ioff = glob.maxblk;
    nsblk = cgetq->cnblk;
  }
  else if (last != 0)
  {
    cgetq = cgetv;
    nsblk = getqq(iblk, &cgetq->cff, cgetq->cwt, ibkptr, ikmin, ivs);
    if (nsblk == 0)
    {
      *idgn = 0;
      return glob.nqn;
    }
    last = 0;
    ioff = 0;
    cgetq->cnblk = nsblk;
    cgetq->cblk = iblk;
    cgetq->csblk = 0;
  }
  else
  {
    cgetq = &cgetv[1];
    ioff = glob.maxblk;
    nsblk = getqq(iblk, &cgetq->cff, cgetq->cwt, &ibkptr[ioff],
                  &ikmin[ioff], &ivs[ioff]);
    if (nsblk == 0)
    {
      *idgn = 0;
      return glob.nqn;
    }
    last = 1;
    cgetq->cnblk = nsblk;
    cgetq->cblk = iblk;
    cgetq->csblk = ioff;
  }
  /*  check for request of size */
  if (indx <= 0)
  {
    *idgn = ibkptr[nsblk + ioff];
    return glob.nqn;
  }
  /*  search for sub-block */
  ix = indx - 1;
  ibgn = cgetq->csblk;
  iend = nsblk + ioff;
  if (ix < ibkptr[ibgn])
    ibgn = ioff;
  for (i = ibgn + 1; i < iend; ++i)
  {
    if (ix < ibkptr[i])
      break;
  }
  cgetq->csblk = ibgn = i - 1;
  /*  assemble quanta */
  ncod = ibkptr[i] - ibkptr[ibgn];
  ldgn = cgetq->cwt[0];
  ngsym = (short)setgsym(cgetq->cwt[3]);
  iqf = cgetq->cff;
  iqsp = getqs(ivs[ibgn], iqf, 0, 0, ixcom, iscom, &ivbase) - 1;
  iv = ixcom[XVIB];
  isym = ixcom[XSYM];
  n = ixcom[XNVAL];
  k = ix - ibkptr[ibgn];
  kq = ikmin[ibgn] + (k << 1);
  if (ngsym >= 3)
  {
    mgsym = (short)MOD(kq - ixcom[XLVAL] + ixcom[XISYM] + ngsym, ngsym);
    if (mgsym != 0)
    {
      if ((mgsym + mgsym) == ngsym)
      { /* B symmetry */
        ldgn = cgetq->cwt[2];
      }
      else
      { /* E symmetry */
        if (glob.esym)
          isym = 4;
        ldgn = cgetq->cwt[1];
      }
    }
  }
  lv = ixcom[XLVAL];
  if (lv != 0)
  { /* l-doubled state */
    if (glob.lsym)
      isym = 4;
    ncod = 1;
    if (kq == 0 && isym == 3)
    {
      if (iv == ivbase)
        ++iv;
      else
        --iv;
    }
  }
  if (glob.nqnn == 2)
  {
    if (ODD(isym))
      kq = -kq;
    ncod = 1;
  }
  else
  {
    if (isym <= 3)
      lower = (n + isym) & 1;
    else /* degenerate state */
      lower = (lv > 0) ? 1 : 0;
    if (glob.oblate)
    {
      iqn[2] = (short)kq;
      kq = n - kq + lower;
      ncod = -ncod;
    }
    else
    {
      iqn[2] = (short)(n - kq + lower);
    }
  }
  iqn[0] = (short)n;
  iqn[1] = (short)kq;
  k = glob.nqnn;
  if (glob.vibfmt)
    iqn[k++] = (short)iv;
  if (cgetq->cwt[4] > maxqn)
  {
    iqn[k++] = (short)iqsp;
    iqn[k] = (short)((iqf + 1) >> 1);
  }
  else
  {
    for (i = 0; k < maxqn; ++i)
    {
      iqn[k++] = (short)((iscom[i] + 1) >> 1);
      if (i == itsym && ngsym > 3)
        i += ngsym - 3;
    }
  }
  ldgn *= iqf + 1;
  if (ldgn < 0)
    ldgn = 0;
  *idgn = ldgn;
  return ncod;
} /* getqn */

/**
 * @brief Calculates direction cosine matrix elements <bra|Phi_Fgz|ket> in a K-basis.
 *        Handles various tensor orders (L=0,1,2) and K changes (DeltaK).
 *        Also incorporates l-doubling and I_tot symmetry effects via operator_flags and selection_mask.
 * @param xbra_qns Quantum numbers for the bra state.
 * @param xket_qns Quantum numbers for the ket state.
 * @param dir_cos_order Tensor order L of the direction cosine operator.
 * @param k_delta_op Overall Delta K for the operator.
 * @param max_elements_out Max number of elements that can be stored in output arrays.
 * @param dir_cos_elements_out Output: Array to store the calculated direction cosine matrix elements (K-dependent part).
 * @param ibra_indices_out Output: Array of bra K-indices (relative to K_begin_bra) for each element.
 * @param iket_indices_out Output: Array of ket K-indices (relative to K_begin_ket) for each element.
 * @param ifup_triangle_flag Controls calculation for upper/lower triangle. If <0 (electric dipole), affects phase.
 * @param operator_flags Bit flags for operator properties (MODD, MLZ from SPAR/SDIP).
 * @param selection_mask 3-bit mask from getmask() indicating allowed K connections (K'=K+DK, K'=DK-K, K'=K-DK).
 * @param is_unit_matrix_type_ptr Output/Input: If input k_avg_val > 0, it's a hint for unit matrix.
 *                                Output is 1 if operator is effectively a unit matrix (K-independent), 0 otherwise.
 * @return int Number of non-zero direction cosine elements calculated. Returns negative of deficit if max_elements_out is too small.
 */
int dircos(xbra, xket, ld, kd, ncmax, direl, ibra, iket, ifup, loff, mask,
           isunit)
    const int *xbra,
    *xket;
const int ld, kd, ncmax, ifup, loff;
int mask;
double *direl;
short *ibra, *iket;
int *isunit;
{ /*  SUBROUTINE TO  CALCULATE DIRECTION COSINE ELEMENTS */
/* XBRA,XKET INTEGER VECTORS CONTAINING: */
/*           0: DIMENSION OF SUB-BLOCK */
/*           1: SYMMETRY CODE FOR BLOCK (0=A,1=Bx,2=By,3=Bz) */
/*           2: 'N' QUANTUM NUMBER */
/*           3: BEGINNING K */
/*           4: L VALUE */
/*           5: VIBRATIONAL QUANTA */
/*           6: SPIN SYMMETRY */
/* L IS TENSOR ORDER OF COSINE */
/* KD IS TOTAL K QUANTIUM CHANGE */
/* DIREL IS A VECTOR OF DIRECTION COSINE ELEMENTS */
/* NCMAX IS MAXIMUM LENGTH OF DIREL */
/* IBRA,IKET INDICES OF MATRIX ELEMENTS IN DIREL (START WITH ZERO) */
/* IFUP =0 IF NO UPPER TRIANGLE TO BE CALC. */
/* (LOFF AND MSYM2) IF OPERATOR HAS ALTERNATE SYMMETRY UNDER ITOT */
/* (LOFF AND MLZ) IF OPERATOR INCLUDES LZ */
/* (LOFF AND MODD) IF OPERATOR IS ODD-ORDER */
/*  L=0  KD=0    DIREL=1. */
/* IFUP <0 IF ELECTRIC DIPOLE */
/* ODD ORDER PARAMETERS ARE IMAGINARY AND */
/* EVEN ARE REAL IF NOT ELECTRIC DIPOLE */
/* FOR ELECTRIC DIPOLE ODD ORDER IS REAL AND EVEN ARE IMAGINARY  */
/*  L=1  KD=0    DIREL=PHI(Z) */
/*  L=1  KD=1    DIREL=PHI(X) OR PHI(Y) */
/*  L=2  KD=0    DIREL=0.5*(3* PHI(Z)**2 -1.) */
/*  L=2  KD=1    DIREL=PHI(Z)*PHI(X)+PHI(X)*PHI(Z)  OR */
/*                      =PHI(Z)*PHI(Y)+PHI(Y)*PHI(Z) */
/*  L=2  KD=2    DIREL=2*(PHI(X)**2 - PHI(Y)**2) OR */
/*                      =PHI(X)*PHI(Y)+PHI(Y)*PHI(X) */
/*     IF KD > L THEN (P-)**(KD-L) OPERATOR MULTIPLIED BY */
/*               L,L OPERATOR */
/*     WANG FUNCTIONS MULTIPLIED BY i**(SYMMETRY) TO MAKE ALL */
/*               ODD ANGULAR MOMENTA IMAG. AND ALL EVEN ONES REAL */
/* IFC = FOURIER SERIES, > 0 IS COS, < 0 IS SIN */
/* MASK = 3-BIT MASK TO ENABLE MATRIX ELEMENT SUB-TYPE */
/*     BEST RESULTS ON SPEED IF N FOR BRA CHANGES SLOWEST */
/***********************************************************************/
#define NFAC_DIRCOS 20
  static double ff[MAXN_DIRCOS], fac[NFAC_DIRCOS];
  double ele, elex, eletmp;
  int kdel, kbra, kbgn, nbgn, kket, isgn, ncase0;
  int ncase, isum, kbra0, kdel2, kket0, kket2, i, k;
  int n, ikdel, kkdel, isbra, nnbra, nsbra, isket, idif;
  int ii, il, kk, lld, ir, ks, nnket, nsket, isrev, nval;
  int iphaz, isym, nbra, nket, kbra0x, kket0x, kavg, kbit;
  BOOL matsym;

  kbra0 = xbra[XKBGN];
  kbra0x = xbra[XIQN] & 1;
  kket0 = xket[XKBGN];
  kket0x = xket[XIQN] & 1;
  nsbra = xbra[XDIM];
  isbra = xbra[XSYM];
  nsket = xket[XDIM];
  isket = xket[XSYM];
  isym = (isbra ^ isket) & 3;
  if ((loff & MODD) != 0)
  {
    if (ncmax == 0)
    { /* initialize */
      fac[0] = 0.5;
      for (lld = 1; lld < NFAC_DIRCOS; ++lld)
      {
        ele = lld - 0.5;
        fac[lld] = fac[lld - 1] * sqrt(lld / ele);
      }
      fac[0] = 1. / fac[1];
      return 0;
    };
    isym = 3 - isym;
  }
  ncase = 0;
  kavg = (*isunit);
  *isunit = 0;
  /****** set up state phases, isum is power of i */
  isrev = isym >> 1;
  /* isym is odd if operator is imaginary */
  isum = ixphase[isket] - ixphase[isbra] + (isym & 1);
  if (ifup < 0) /* correct to make dipoles real */
    ++isum;
  if (ODD(isum) && xbra[XVIB] < xket[XVIB])
    isum += 2; /* imaginary vibration phase assumed */
  isgn = isum >> 1;
  if ((loff & MLZ) != 0 && xket[XLVAL] < 0)
    ++isgn; /* Lz operator found */
  /*********** TRY SCALAR ******************************/
  if (kd == 0 && ld == 0 && (mask & 5) != 0)
  {
    mask &= 2;
    ii = 0;
    ncase = nsket;
    n = nsbra;
    kdel = kbra0 - kket0;
    kk = kbra0;
    if (kdel != 0)
    {
      if (kdel > 0)
      {
        kdel = kdel >> 1;
        ncase -= kdel;
        if (ncase < 0)
          ncase = 0;
        kk = kbra0;
      }
      else
      {
        ii = (-kdel) >> 1;
        kdel = -ii;
        n -= ii;
        if (n < 0)
          n = 0;
        kk = kket0;
      }
    }
    if (ncase > n)
      ncase = n;
    if (ncase > 0)
    {
      idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      if (kavg > 0)
      { /* single K value */
        k = kavg - kk;
        if (k < 0 || ODD(k))
          return 0;
        k = k >> 1;
        if (k >= ncase)
          return 0;
        ii += k;
        kk += k;
        ncase = 1;
      }
      if (kk != 0)
      { /* kk = max beginning k */
        mask = 0;
      }
      else if (kbra0x == kket0x)
      {
        if (kbra0x == 0)
          mask = 0;
        kk = 2;
      }
      else
      { /* kk = 0 */
        mask = 0;
      }
      ele = 1.;
      if (ODD(isgn + isrev))
        ele = -ele;
      *isunit = 1;
      for (k = 0; k < ncase; ++k)
      {
        direl[k] = ele;
        ibra[k] = (short)ii;
        iket[k] = (short)(ii + kdel);
        ++ii;
      }
      if (kk == 0)
      {
        direl[0] *= fac[0];
        *isunit = 0;
      }
    }
    if (mask == 0)
      return ncase; /* no quasi-diagonal delta K = 0 */
  } /************** end Scalar *********/
  nbra = xbra[XNVAL];
  nket = xket[XNVAL];
  kdel = kd;
  ikdel = kdel - ld;
  if (ikdel <= 0)
  {
    kkdel = kdel;
  }
  else
  {
    k = kbra0 + ikdel + ((nsbra - 1) << 1);
    ffcal(nbra, k, ff);
    kkdel = ld;
  }
  kdel2 = kkdel + kkdel;
  lld = ld + ld;
  nnbra = nbra + nbra;
  nnket = nket + nket;
  if (ld > 0)
  {
    ele = (double)(nnbra + 1);
    if (nnket != nnbra)
      ele = sqrt(ele * (nnket + 1));
    if (kkdel != 0)
      ele *= fac[ld];
    isgn += nbra + kket0 + kkdel;
  }
  else
  {
    ele = 1.;
    if (kdel != 0)
      ele = 0.5;
  }
  matsym = TRUE;
  k = mask & 5;
  if (kdel == 0)
  {
    if (k != 0)
      mask = (mask & 2) + 1;
    matsym = FALSE;
  }
  else
  {
    if (k == 1 || k == 4)
      matsym = FALSE;
  }
  if (matsym)
  {
    for (k = 1; k <= 5; ++k)
    {
      if (xbra[k] != xket[k])
      {
        matsym = FALSE;
        break;
      }
    }
  }
  if (ODD(isgn + isrev))
    ele = -ele;
  kbit = kavg + kavg - (kdel & 1);
  ncase0 = 0;
  for (iphaz = 1; iphaz <= 4; iphaz = iphaz << 1)
  {
    /* iphaz = 1,2,4 */
    if (mask < iphaz)
      break;
    if (iphaz != 2)
    {
      /* iphaz = 1: +KDEL, iphaz = 4: -KDEL */
      /*  START CALCULATION OF  <K ....K-KDEL> */
      /*     KBRA=KBRA0+ (IL+ i)*KINC   IL >= 0 */
      /*     KKET=KKET0+ (IR+ i)*KINC   IR >= 0, i=0,... */
      /*     KBRA=KKET+KDEL */
      if ((mask & iphaz) == 0)
        continue;
      if (iphaz == 4)
      {
        if (matsym)
          break;
        kdel = -kdel;
        kkdel = -kkdel;
      }
      isum = kdel + kket0 - kbra0; /* ISUM = (IL-IR)*KINC */
      if (isum >= 0)
      {
        il = isum >> 1;
        ir = 0;
        kbra = kdel + kket0;
        kket = kket0;
      }
      else
      {
        il = 0;
        ir = (-isum) >> 1;
        kbra = kbra0;
        kket = kbra0 - kdel;
      }
      n = nsket - ir;
      i = nsbra - il;
      if (n > i)
        n = i;
      if (n <= 0)
        continue;
      if (kavg > 0)
      {
        kk = (kbit - kbra - kket) >> 1;
        if (kk < 0 || ODD(kk))
          continue;
        kket += kk;
        kbra += kk;
        kk = kk >> 1;
        if (kk >= n)
          continue;
        il += kk;
        ir += kk;
        n = 1;
      }
      nbgn = ncase;
      ncase += n;
      idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      kket2 = kket + kket;
      ncase0 = n;
      if (ld == 0)
      {
        for (n = nbgn; n < ncase; ++n)
        {
          direl[n] = ele;
          ibra[n] = (short)il;
          iket[n] = (short)ir;
          ++il;
          ++ir;
        }
      }
      else
      {
        for (n = nbgn; n < ncase; ++n)
        {
          kk = -kdel2 - kket2;
          direl[n] = ele * c3jj(nnbra, lld, nnket, kk, kdel2, kket2);
          kket2 += 4;
          ibra[n] = (short)il;
          iket[n] = (short)ir;
          ++il;
          ++ir;
        }
      }
      il = ir = 1;
      if (kbra == 0)
        il = kbra0x;
      if (kket == 0)
        ir = kket0x;
      if (ODD(il + ir)) /*  CORRECT FOR K=0 */
        direl[nbgn] *= fac[0];
      if (ikdel > 0)
      { /* ADD ON P- OR P+ OPERATOR TO LEFT */
        kbgn = kbra;
        if (kdel > 0)
          kbgn -= ikdel;
        for (n = nbgn; n < ncase; ++n)
        {
          eletmp = direl[n];
          kk = kbgn + ikdel;
          for (k = kbgn; k < kk; ++k)
          {
            eletmp *= ff[k];
          }
          direl[n] = eletmp;
          kbgn += 2;
        }
      }
    }
    else
    { /* iphaz == 2 */
      /*  DO QUASI-DIAGONAL PART */
      /*   <KDEL-K ..... K>   K=1,...,KDEL-1 */
      /*       KBRA=KBRA0+ (IL - i)*KINC  IL <= n-1 */
      /*       KKET=KKET0+ (IR + i)*KINC  IR >= 0  , i=0 .. n-1 */
      /*       KDEL=KKET+KBRA */
      kdel2 = -kdel2;
      if (ODD(isrev))
        ele = -ele;
      if ((mask & 2) == 0)
        continue;
      isum = kdel - kbra0 - kket0;
      if (isum < 0)
        continue;
      isum = isum >> 1; /* isum = (IL + IR) */
      ir = 0;           /* IR = 0 */
      if (isum >= nsbra)
        ir = isum - nsbra + 1;
      else if (kket0x == 0)
        ++ir;
      n = isum + 1;
      if (n > nsket)
        n = nsket;
      else if (kbra0x == 0)
        --n;
      n -= ir;
      if (n <= 0)
        continue;
      il = isum - ir;
      if (ifup == 0)
        n = (il - ir + 2) >> 1;
      kket = kket0 + (ir << 1);
      kbra = kdel - kket;
      nbgn = ncase;
      ncase += n;
      idif = ncmax - ncase;
      if (idif < 0)
        return idif;
      if (ODD2(isbra + nnbra))
      {
        elex = -ele;
      }
      else
      {
        elex = ele;
      }
      n = nbgn;
      isgn = loff & MFCM;
      for (nval = nbgn; nval < ncase; ++nval)
      {
        kk = 0;
        if (kavg > 0)
        {
          kk = kket - kbra;
          if (kk < 0)
            kk = -kk;
          kk -= kbit;
        }
        else if (kbra == kket)
        {
          kk = isgn; /* skip if sin(0) */
        }
        if (kk == 0)
        {
          eletmp = elex;
          if (isgn != 0 && kbra > kket)
            eletmp = -eletmp;
          if (ld != 0)
          {
            kket2 = kket + kket;
            ks = -kdel2 - kket2;
            eletmp *= c3jj(nnbra, lld, nnket, ks, kdel2, kket2);
          }
          for (kk = 0; kk < ikdel; ++kk)
          {
            k = kbra - kk;
            if (k <= 0)
            {
              eletmp *= ff[-k];
            }
            else
            {
              eletmp *= ff[k - 1];
            }
          }
          direl[n] = eletmp;
          ibra[n] = (short)il;
          iket[n] = (short)ir;
          ++n;
        }
        --il;
        ++ir;
        kbra -= 2;
        kket += 2;
      } /* loop over elements */
      ncase = n;
    }
  } /* loop over 3 cases */
  if (matsym && ncase0 > 0 && ifup != 0)
  {
    k = ncase;
    ncase += ncase0;
    idif = ncmax - ncase;
    if (idif < 0)
      return idif;
    for (i = 0; i < ncase0; ++i, ++k)
    {
      direl[k] = direl[i];
      ibra[k] = iket[i];
      iket[k] = ibra[i];
    }
  }
  return ncase;
} /* dircos */

/**
 * @brief Calculates factors related to N_K = sqrt(N(N+1)-K(K-1)) for P+- operators.
 *        Stores sqrt( N(N+1) - K(K-1) ) for K from 0 up to k_val_limit.
 *        Used by dircos when operator involves P+ or P- type terms (DK_op > L_op).
 * @param n_val_for_factor N quantum number.
 * @param k_val_limit Maximum K value for which to calculate the factor.
 * @param output_factors_array Output: Array to store the calculated factors ff[K] = sqrt(N(N+1)-K(K-1)).
 * @return int Always 0.
 */
int ffcal(nff, kff, ff)
    const int nff,
    kff;
double *ff;
{
  static int nlast = -1;
  static int klast = 0;
  static double sq;
  int k;

  if (nff != nlast)
  {
    nlast = nff;
    sq = (double)nff;
    sq *= nff + 1;
    k = 0;
  }
  else
  {
    if (klast >= kff)
      return 0;
    k = klast;
  }
  klast = kff;
  if (klast > nff)
    klast = nff;
  while (k < klast)
  {
    ff[k] = sqrt(sq);
    ++k;
    sq -= k + k;
  }
  return 0;
} /* ffcal */
