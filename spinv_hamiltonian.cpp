/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 3 March 2004 */
/*   28 Sept.  2005: fix phase bug in getzitot for reduced matrix elements */

/* split from original spinv.c */
/* Functions related to the construction and manipulation of the Hamiltonian matrix */

/*
hamx(): The main Hamiltonian calculation function.
specop(): Applies Euler series operator transformations.
specfc(): Applies Fourier coefficient scaling.
sznzfix(), sznzop(): Handle S<sub>z</sub>N<sub>z</sub> and related N-changing operators.
rmatrx(): N-dependent corrections for reduced matrix elements.
symnsq(): Applies N(N+1) and N.S operator scaling.
symksq(): Applies K<sup>2</sup> operator scaling.
dpmake(): Calculates derivative contributions (d<H>/dP).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "calpgm.h"
#include "spinit.h"
#include "spinv_internal.h"
#include "SpinvContext.hpp"
#include "CalError.hpp"

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
int hamx(struct SpinvContext *ctx, const int iblk, const int nsize, const int npar,
         const bcd_t *idpar, const double *par, double *egy, double *t, double *dedp,
         double *pmix, BOOL ifdump)
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
  int ipar, ncos, ispn, jspn, nsqj, iiwt[5], ndmd, alpha, mkd, njqt; /* ipar: parameter index; ncos: number of non-zero direction cosine elements; ispn,jspn: packed spin state for bra/ket; nsqj: power of N(N+1); iiwt: cached weights for current block; ndmd: nsize+1 for diagonal indexing; alpha: Itot symmetry component; mkd: mask for dircos */
  int i, ii, lt, n_sub_block_size, ibase, jbase, kbgni, kbgnj, nsblk, idflags, ipbase; /* i,ii: loop counters; lt: tensor L order; n_sub_block_size: current sub-block size; ibase,jbase: base index for bra/ket sub-block; kbgni,kbgnj: K_begin for bra/ket; nsblk: number of sub-blocks in current F block; idflags: flags from SPAR struct; ipbase: base parameter index for constrained params */
  int kd, ld, nd, ni, nj, ivbase, jvbase, ivmin, ifc, iz, jz, npair, ldel; /* kd: Delta K; ld: dir.cos. L order; nd: N_bra - N_ket; ni,nj: N_bra, N_ket; ivbase,jvbase: vib. index for bra/ket; ivmin: min(ivbase,jvbase); ifc: Fourier coeff index; iz,jz: matrix indices; npair: number of spin couplings; ldel: Delta l */
  int si1, si2, iff, ijd, ikq, ins, neuler, sznz, ixx, jxx, kl, isgn, isunit; /* si1,si2: spin indices; iff: 2*F; ijd: ixx-jxx; ikq: K^2 power; ins: N.S power; neuler: Euler type; sznz: SzNz type; ixx,jxx: sub-block loop indices; kl: loff flags for dircos; isgn: sign factor for matrix element; isunit: flag if operator is unit matrix type */
  int kavg, nqnsep; /* kavg: K_average for specific K parameters; nqnsep: result from hdiag, number of separated states or error */
  unsigned int ivsym, lastvv, ivcmp; /* Packed vibrational and symmetry identifiers for fast parameter lookup */
  BOOL oldpar, isegy0, roll, parskp, isneg, first; /* oldpar: if previous param shared cosine part; isegy0: if current operator contributes to E0; roll: if K-roll detected; parskp: skip current parameter; isneg: if oblate hamiltonian needs negation; first: first pass (calculating H) vs second pass (derivatives) */
  BOOL newblk, firstpar; /* newblk: if new vibrational pair, re-eval Euler denoms; firstpar: first parameter for a given vib pair and symmetry */

  if (ctx->ndmx <= 0) { /* Check if work arrays are allocated */
    throw NumericError("working vectors not allocated", CalErrorCode::WorkingVectorTooShort);
  }
  dbar = &ctx->wk[ctx->ndmx]; /* dbar points to a section of the work array wk, used for accumulating derivatives for a single operator */
  ndm = nsize;      /* Store nsize in ndm (leading dimension for matrix t) */
  ndmd = nsize + 1; /* Used for indexing diagonal elements t[i+i*ndm] as t[i*ndmd] */
  isneg = (ctx->glob.oblate && !ifdump); /* Hamiltonian needs to be negated for oblate rotors if not just dumping */
  ctx->cgetv[0].cblk = 0;                /* Invalidate cache for getqn for state 0 */
  ctx->cgetv[1].cblk = 0;                /* Invalidate cache for getqn for state 1 */
  /*     get F and sub-block structure */
  itau = &(ctx->ikmin)[ctx->glob.maxblk]; /* itau points to an area after ikmin, used to store sorting keys (tau or K values) */
  nsblk = getqq(ctx, iblk, &iff, iiwt, ctx->ibkptr, ctx->ikmin, ctx->ivs); /* Get sub-block structure for current F block (iblk) */
  /* zero hamiltonian matrix t */
  dclr(nsize, nsize, t, 1);
  /* zero derivative matrix dedp */
  n_sub_block_size = ctx->glob.nfit; /* Number of parameters to fit */
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
      ibase = ctx->ibkptr[ixx]; /* Starting row/col index for this 'bra' sub-block */
      n_sub_block_size = ctx->ibkptr[ixx + 1] - ibase; /* Size of this 'bra' sub-block */
      kbgni = ctx->ikmin[ixx]; /* Starting K for this 'bra' sub-block */
      ispn = ctx->ivs[ixx];    /* Packed (Vib,Sym,SpinPattern) for 'bra' */
      getqs(ctx, ispn, iff, n_sub_block_size, kbgni, ctx->ixcom, ctx->iscom, &ivbase); /* Get quantum numbers for 'bra' state */
      ni = ctx->ixcom[XNVAL]; /* N quantum number for 'bra' */
      /* set up itau for sorting (tau = Ka-Kc like proxy, or K for projection) */
      kd = (ni + ctx->ixcom[XSYM] + 1) & 1; /* Parity component for tau */
      kd += kbgni + kbgni;             /* K component for tau */
      kd = kd + kd;                    /* Scale it up */
      if (ivbase != ctx->ixcom[XVIB]) /* Check if it's the upper state of an l-doublet pair based on how getqs returns XVIB relative to the actual state index ivbase */
        ++kd; /* Adjust tau for upper l-doubled state */
      itau[ixx] = (short) kd; /* Store sorting key for this sub-block */

      if (ni != nold) {  /* If N changed, re-calculate powers of N(N+1) */
        nold = ni;
        dn = (double) ni;
        sqnn = dn * (ni + 1); /* N(N+1) */
        sqj[0] = 1.0; /* Zeroth power */
        for (i = 1; i <= ctx->nsqmax; ++i)
        { /* Calculate higher powers */
          sqj[i] = sqj[i - 1] * sqnn;
        }
      }

      for (jxx = 0; jxx <= ixx; ++jxx) { /* Loop over 'ket' sub-blocks (only jxx <= ixx for lower triangle + diagonal) */
        ijd = ixx - jxx; /* Difference in sub-block indices */
        jbase = ctx->ibkptr[jxx]; /* Starting row/col index for this 'ket' sub-block */
        n_sub_block_size = ctx->ibkptr[jxx + 1] - jbase; /* Size of this 'ket' sub-block */
        kbgnj = ctx->ikmin[jxx];                         /* Starting K for this 'ket' sub-block */
        jspn = ctx->ivs[jxx];    /* Packed (Vib,Sym,SpinPattern) for 'ket' */
        getqs(ctx, jspn, iff, n_sub_block_size, kbgnj, ctx->jxcom, ctx->jscom, &jvbase); /* Get quantum numbers for 'ket' state */
        nj = ctx->jxcom[XNVAL];                                                     /* N quantum number for 'ket' */
        nd = ni - nj;      /* Delta N */
        if (nd > ctx->ndmax && ctx->ndmax >= 0)
          continue; /* Skip if Delta N too large (ndmax < 0 allows all Delta N) */

        ivmin = (ivbase < jvbase) ? ivbase : jvbase; /* Minimum of the two vibrational indices */
        /* Create a packed identifier for the vibrational pair and overall block symmetry */
        ivsym = (unsigned int) (ivbase + jvbase + ivmin * ctx->glob.vibfac) << 2;
        ivcmp = ivsym + 3; /* Upper bound for checking Euler parameters which might be only vv' specific */
        ivsym += blksym(ctx->ixcom, ctx->jxcom); /* Add block symmetry (A,Bx,By,Bz) */

        if (lastvv != ivcmp) { /* If vibrational pair changed, Euler denominators might need re-init */
          newblk = TRUE; lastvv = ivcmp;
        }
        isegy0 = FALSE; firstpar = TRUE; /* isegy0: current op contributes to E0; firstpar: first param for this (v,v',sym) */

        for (spar_now = ctx->spar_head[ivmin]; TRUE; spar_now = spar_now->next)
        { /* Loop over parameters relevant to this vibrational pair */
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
                  iz = ctx->idx[0]; /* Only one set of indices for unit matrix */
                  jz = ctx->jdx[0];
                  pt = &t[jz * ndm]; /* Pointer to start of column jz */
                  for (i = 0; i < ncos; ++i) { /* ncos is actual # elements for unit matrix */
                    pt[iz] += pbar; /* Add to H[iz, jz] (actually H[row_idx[i], col_idx[i]])*/
                    iz += ndmd; /* Move to next element assuming it's effectively on a diagonal of a submatrix if isunit represents block diagonal structure */
                                /* This indexing for unit matrix seems specific; idx[0], jdx[0] are start, ncos is length of this unit block along diagonal */
                  }
                } else { /* Operator has K-dependent part stored in wk array */
                  for (i = 0; i < ncos; ++i) {
                    iz = ctx->idx[i]; /* Row index for this element */
                    jz = ctx->jdx[i]; /* Column index for this element */
                    t[iz + jz * ndm] += pbar * ctx->wk[i]; /* Add PBAR * K_dependent_factor */
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
              sznzfix(ctx, sznz, ni, nj, ctx->ixcom, ctx->jxcom, ctx->iscom, ctx->jscom); /* Reset N values in ixcom/jxcom if changed by SzNz */
            /* Parse the BCD parameter ID into its components */
            kl = idpars(spar_now, &ikq, &neuler, &lt, &ld, &kd, &ins,
                        &si1, &si2, &sznz, &ifc, &alpha, &ldel, &kavg, &njqt);
            if (sznz > 0) { /* If this is an SzNz type operator */
              sznz = sznzfix(ctx, sznz, ni, nj, ctx->ixcom, ctx->jxcom, ctx->iscom, ctx->jscom); /* Modify N in ixcom/jxcom, returns new sznz state */
              if (sznz == 0) continue;    /* Operator invalid for these N values, get next parameter */
            }
            if (kavg > 0 && kavg > ni && kavg > nj) continue; /* If specific K_avg is outside range of N_bra, N_ket */

            mkd = getmask(ctx, ctx->ixcom, ctx->jxcom, kd, ldel, kl, alpha); /* Get selection mask for dircos based on symmetry */
            if (mkd == 0) continue; /* Operator forbidden by symmetry/selection */

            npair = getll(ctx, 0, ld, lt, kd, si1, si2, ctx->lscom, ctx->iscom, ctx->jscom); /* Get spin coupling tensor orders */
            if (npair < 0) continue; /* Tensor coupling rules violated */

            if (ncos < 0) { /* If K-dependent part needs recalculation */
              isunit = kavg; /* kavg might indicate a unit matrix type if > 0 and nofc=1 */
              /* Calculate direction cosine matrix elements (K-dependent part) */
              ncos = dircos(ctx, ctx->ixcom, ctx->jxcom, ld, kd, ctx->ndmx, ctx->wk, ctx->idx,
                            ctx->jdx, ijd, kl, mkd, &isunit);
              if (ncos <= 0) {
                if (ncos == 0) continue; /* No non-zero elements, get next parameter */
                throw NumericError("DIRCOS WORKING VECTOR TOO SHORT IN HAMX", CalErrorCode::WorkingVectorTooShort);
              }
              if (TEST(idflags, MNOUNIT)) /* If flag says it's not a unit matrix, ensure isunit is 0 */
                isunit = 0;
              isgn = 1;
              if (isunit != 0 && ctx->wk[0] < 0.) /* For unit matrix type, if scaling factor is negative */
                isgn = -1;

              if (neuler != 0) { /* If Euler series operator */
                /* Apply Euler series transformation to ctx->wk, idx, jdx */
                /* par[ipar] is the Euler denominator 'a' or 'b' for this specific neuler type */
                ncos = specop(neuler, &newblk, &nsqj, &ikq, kbgni, kbgnj,
                              ni, nj, ncos, ctx->wk, ctx->idx, ctx->jdx, par[ipar]);
                if (ncos == 0) continue; /* Euler op became zero */
              }
              if (ifc != 0) /* If Fourier coefficient operator */
                specfc(ctx, ifc, ivbase, jvbase, kd, kbgni, kbgnj, ncos,
                       ctx->wk, ctx->idx, ctx->jdx); /* Apply Fourier scaling to wk */
              if (sznz != 0) /* If SzNz type operator */
                sznzop(ctx, ni, nj, kbgni, kbgnj, ctx->iscom, ctx->jscom, ncos, ctx->wk, ctx->idx, ctx->jdx); /* Apply SzNz scaling to wk */
              if (ikq > 0) { /* If K^2 dependent operator */
                ncos = symksq(ikq, kbgni, kbgnj, ncos, ctx->wk, ctx->idx, ctx->jdx); /* Apply K^2 scaling */
                if (ncos == 0) continue;
              }
              /* Adjust indices if not in the first sub-block of the matrix */
              if (ibase != 0) {
                for (i = 0; i < ncos; ++i) {
                  ctx->idx[i] = (short)(ctx->idx[i] + ibase);
                }
              }
              if (jbase != 0) {
                for (i = 0; i < ncos; ++i) {
                  ctx->jdx[i] = (short)(ctx->jdx[i] + jbase);
                }
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
                dpmake(nsize, dbar, t, ncos, ctx->wk, ctx->idx, ctx->jdx, isunit);
              }
            } /* End K-dependent part calculation (if ncos < 0) */

            if (si2 < 0 && nd == 0 && ni ==0) continue;  /* Avoids issues with commutator with N*N when N=0 */

            /* Calculate K-independent part of operator matrix element */
            zpar = rmatrx(ld, lt, ctx->ixcom, ctx->jxcom); /* N-dependent corrections for reduced matrix elements */
            tensor(ctx, &zpar, ctx->iscom, ctx->jscom, ctx->lscom, ctx->ismap, npair, alpha); /* Spin-dependent part (Clebsch-Gordan etc.) */

            if (sznz < 0) /* If SzNz type needs finalization from idpars */
              sznz = sznzfix(ctx, sznz, ni, nj, ctx->ixcom, ctx->jxcom, ctx->iscom, ctx->jscom); /* Reset N in ixcom/jxcom if changed */

            symnsq(ctx, nsqj, ins, ctx->iscom, ctx->jscom, &zpar); /* N(N+1) and N.S dependence */
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

          i = ctx->ipder[ipar]; /* Get derivative index for this parameter */
          ipbase = ipar;
          if (i < 0) { /* If parameter is constrained (i < 0 points to main parameter) */
            ipbase = -1 - i; /* Index of the parameter it's constrained to */
            i = ctx->ipder[ipbase]; /* Derivative index of the main parameter */
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
      if (ctx->glob.idiag < 0 || nsize <= 1) { /* If no diagonalization requested or trivial 1x1 block */
        nd = nsize + 1;
        dcopy(nsize, t, nd, egy, 1); /* Copy diagonal elements of H to egy */
        for (i = 0; i < nsize; ++i) { /* Create identity eigenvector matrix */
          pt = &t[i * ndm];
          memset(pt, 0, sizeof(double) * nsize); /* Zero out column */
          pmix[i] = pt[i] = 1.;          /* Eigenvector is unit vector, pmix=1 */
        }
        if (ctx->glob.idiag == 4 && nsize == 1) /* Special case for K-sort of 1x1 */
          pmix[0] = (double) kbgni; /* pmix stores K value */
      } else { /* Diagonalization is requested and nsize > 1 */
        if (ctx->glob.idiag == 0) { /* Standard energy sort, check for K-roll before diagonalization */
          roll = kroll(nsize, t, nsblk, ctx->ibkptr, ctx->ikmin);  /* Modifies t if roll-over found */
        } else if (ctx->glob.idiag == 2 || ctx->glob.idiag == 5) { /* Sort by input H diagonal order */
          nd = nsize + 1;
          dcopy (nsize, t, nd, dbar, 1); /* Save original diagonal of H in dbar for later sorting */
        }
        /* diagonalize Hamiltonian matrix t, eigenvalues in egy, eigenvectors in t (overwriting H) */
#if HAM_DEBUG /* Debug block to check matrix before/after diagonalization */
        nd = nsize * nsize;
        pscr = (double *) calalloc((size_t)nd * sizeof(double)); /* Scratch space */
        if (nd > 0) pscr[0] = t[0]; /* To use pscr if nd > 0 */
        dcopy(nd, t, 1, pscr, 1); /* Save original H */
#endif
        nqnsep = hdiag(nsize, nsize, t, egy, pmix, ctx->iqnsep); /* Diagonalize */
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
          throw NumericError("diagonalization failure in hamx", CalErrorCode::DiagonalizationFailed);
        }
      	switch (ctx->glob.idiag) { /* Sorting options for eigenvalues/vectors */
      	default:  /* idiag = 0: energy sort of subblock (Wang sub-blocks are sorted internally) */
          ordblk(nsize, nsize, ctx->iqnsep, t, egy, ctx->ibkptr, pmix, ctx->idx); /* idx is scratch */
          break;
        case 1:  /* projection sort of full block (sort by <K> or similar from pmix) */
          for (i = 0; i <= nsize; ++i) { /* Create identity permutation for full block sort */
            ctx->jdx[i] = (short)i;
          }
          ordblk(nsize, nsize, ctx->iqnsep, t, egy, ctx->jdx, pmix, ctx->idx); /* idx is scratch */
          break;
        case 2: /* Energy sort within Wang sub-blocks, but ensure it follows original H diagonal order if energies are close */
          i = ordblk(nsize, nsize, ctx->iqnsep, t, egy, ctx->ibkptr, pmix, ctx->idx); /* Initial sort by energy within sub-blocks */
          if (ODD2(i)) { /* ODD2(i) might mean some degeneracy or near-degeneracy was handled by ordblk */
            ordham(nsize, ctx->iqnsep, dbar, ctx->ibkptr, ctx->jdx); /* Get permutation jdx based on original H diagonal (dbar) */
            fixham(nsize, nsize, t, egy, pmix, ctx->jdx);       /* Apply this permutation */
          }
          break;
        case 3: /* Sort by tau = Ka-Kc like proxy, within vibrational and spin sub-blocks */
        case 4: /* Sort by K (projection of K^2), within vibrational and spin sub-blocks */
        case 5: /* Sort by original H diagonal, within vibrational and spin sub-blocks */
          vbkptr = &(ctx->ibkptr)[ctx->glob.maxblk]; /* Use area after ibkptr for vibrational sub-block pointers */
          nd = 0;
          ivbase = -1; /* Current vibrational state index being processed */
          for (i = 0; i < nsblk; ++i) { /* Determine boundaries of vibrational sub-blocks */
            ixx = ctx->ibkptr[i]; /* Start of current (Wang*spin) sub-block */
            /* Group sub-blocks by common vibrational state (and l-doublet component) */
            jvbase = (int) ((unsigned) ctx->ivs[i] >> 2) - (itau[i] & 1); /* Effective vib index for grouping l-doublets */
            if (ivbase != jvbase) { /* If new vibrational group starts */
              vbkptr[nd++] = (short) ixx; /* Store start index of this new group */
              ivbase = jvbase;
            }
          }
          vbkptr[nd] = (short) nsize; /* Mark end of last vibrational group */
          i = ordblk(nsize, nsize, ctx->iqnsep, t, egy, vbkptr, pmix, ctx->idx); /* Initial energy sort within these vib groups */
          if ((i & 2) == 0 /* TODO INCLUDE?: && glob.idiag != 4 */) /* If no issues from ordblk, and not K-sort, we might be done */
            break;                            /* For idiag=4 (bestk), always proceed */
          if (ctx->glob.idiag == 5){ /* Sort by original H diagonal within vib groups */
            ordham(nsize, ctx->iqnsep, dbar, vbkptr, ctx->idx); /* Get permutation based on H_diag (dbar) within vib groups */
            fixham(nsize, nsize, t, egy, pmix, ctx->idx);       /* Apply permutation */
            break;
          }
          /* For idiag=3 (tau) or idiag=4 (K) */
          for (i = 0; i < nsblk; ++i) { /* Fill idx with sorting keys (tau values stored in itau earlier) */
            ixx = itau[i];
            iz = ctx->ibkptr[i + 1];
            for (ii = ctx->ibkptr[i]; ii < iz; ++ii) {
              ctx->idx[ii] = (short)ixx; /* Assign sorting key for each basis state */
              ixx += 8; /* Increment key for next basis state in sub-block, arbitrary large step */
            }
          }
          if (ctx->glob.idiag == 4) { /* Sort by K values */
            /* bestk calculates K expectation values into wk, gets permutation jdx, applies to t, egy, pmix */
            bestk(nsize, nsize, ctx->iqnsep, vbkptr, ctx->idx, ctx->jdx, t, egy, pmix, ctx->wk);
            dcopy(nsize, ctx->wk, 1, pmix, 1);       /* Store K expectation values in pmix */
            fixham(nsize, nsize, t, egy, pmix, ctx->jdx); /* Apply final permutation for K sort */
            break;
          }
          /* idiag == 3 (tau sort) */
          for (i = 0; i < nsize; ++i) { /* Copy tau keys from idx to wk (double) */
            ctx->wk[i] = (double)ctx->idx[i];
          }
          ordham(nsize, ctx->iqnsep, ctx->wk, vbkptr, ctx->idx); /* Get permutation based on tau values (in wk) */
          fixham(nsize, nsize, t, egy, pmix, ctx->idx);     /* Apply permutation */
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
      	memset(egy, 0, sizeof(double) * nsize); /* Zero out energy array */
	      for (ipar = npar - 1; ipar >= 0; --ipar) { /* Loop over all parameters */
          i = ctx->ipder[ipar];                    /* Get derivative index */
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
int specop(const int neuler, BOOL *newblk, int *nsqj, int *ikq,
           const int ksi, const int ksj, const int ni, const int nj, const int ncos,
           double *wk, const short *ix, const short *jx, const double par)
{
#define NSPOP 5 /* Number of Euler series denominators (a_n, b_n up to n=4, since ist = (neuler-2)/2) */
  static double aden[NSPOP], bden[NSPOP]; /* Denominators 'a' and 'b' for Euler series */
  double akisq, akjsq, anisq, anjsq, da, dd, di, dj, xi, xj; /* K^2_bra, K^2_ket, (N(N+1)-K^2)_bra, (N(N+1)-K^2)_ket, a-d, d, denominator_bra, denominator_ket, element_bra_part, element_ket_part */
  double sqi, sqj, di0, dj0; /* N(N+1)_bra, N(N+1)_ket, (1+d*N(N+1))_bra, (1+d*N(N+1))_ket */
  int i_elem, nkq, nkq0, nnq, nnq0, ki, kj, ist; /* i_elem: loop counter; nkq,nkq0: K^2 power current,original; nnq,nnq0: N(N+1) power current,original; ki,kj: K_bra,K_ket for current element; ist: Euler series index (0 to NSPOP-1) */
  unsigned int ipwr; /* Loop counter for applying powers */

  if (*newblk)
  { /* If new vibrational block, reset denominators */
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
int specfc(struct SpinvContext *ctx, const int ifc, const int iv, const int jv,
           const int kdel, const int ksi, const int ksj, const int ncos,
           double *wk, const short *ix, const short *jx) // TODO: replace ctx parameter with glob.g12
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

  is_g12_term_active = ctx->glob.g12 & fourier_order_n; /* Check if G12 symmetry term (-1)^K is active for this Fourier order */
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
int sznzfix(struct SpinvContext *ctx, const int sznz, const int ni, const int nj,
            int *ixcom, int *jxcom, int *iscom, int *jscom) // TODO remove ctx param with nspin
{
  if (sznz <= 0) { /* Negative or zero sznz means reset N to original values */
    if (sznz <= -4) { /* Reset ket side (nj) */
      jscom[ctx->nspin] = nj << 1; /* jscom[nspin] stores 2*N_ket */
      jxcom[XNVAL] = nj;
    } else { /* Reset bra side (ni) or both if sznz == 0 (no-op but resets state) */
      iscom[ctx->nspin] = ni << 1; /* iscom[nspin] stores 2*N_bra */
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
      iscom[ctx->nspin] -= 2; /* iscom[nspin] is 2*N_bra, so decrement by 2*(1) */
      --ixcom[XNVAL];    /* Decrement N_bra in quantum block */
      break;
    case 3: /* Operator like N_+, changes N_bra to N_bra+1 */
      if (ni < 0) /* N_bra cannot be negative */ /* Original check ni <=0, but N=0 can be raised */
        return 0;
      iscom[ctx->nspin] += 2;
      ++ixcom[XNVAL];
      break;
    case 4: /* Operator like N_-, changes N_ket to N_ket-1 (from ket perspective) */
      if (nj <= 0) /* N_ket cannot be 0 or less to decrease */ /* Original check nj <=1 */
        return 0;
      jscom[ctx->nspin] -= 2;
      --jxcom[XNVAL];
      break;
    case 5: /* Operator like N_+, changes N_ket to N_ket+1 (from ket perspective) */
      if (nj < 0) /* N_ket cannot be negative */ /* Original check nj <= 0 */
        return 0;
      jscom[ctx->nspin] += 2;
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
int sznzop(struct SpinvContext *ctx, const int ni, const int nj, const int ksi,
           const int ksj, const int *iscom, const int *jscom, const int ncos,
           double *wk, const short *ix, const short *jx) // TODO: replace ctx parameter with nspin
{
  int iflg_n_change_type, i_elem, k_val, kk_val, n_bra_times_2, n_ket_times_2, n_bra_orig_times_2, n_ket_orig_times_2, j_qn_times_2, spin_qn_times_2; /* Renamed iflg, i, k, kk, nni, nnj, nni0, nnj0, jj, iis */
  double val_bra, val_ket, current_wk_element, sum_contrib; /* Renamed vali, valj, twk, sum */

  if (ncos <= 0) /* No elements to operate on */
    return 0;

  n_bra_times_2 = ni << 1; /* Current 2*N_bra (possibly modified by sznzfix) */
  n_bra_orig_times_2 = iscom[ctx->nspin];      /* Original 2*N_bra from spin table */
  val_bra = (double) (n_bra_orig_times_2 + 1); /* Factor related to (2N_orig+1) */

  n_ket_times_2 = nj << 1; /* Current 2*N_ket */
  n_ket_orig_times_2 = jscom[ctx->nspin];      /* Original 2*N_ket */
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
    spin_qn_times_2 = jscom[ctx->nspin + 1]; /* 2*S_ket (first spin) */
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
    spin_qn_times_2 = iscom[ctx->nspin + 1]; /* 2*S_bra */
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
 * @brief Calculates N-dependent correction factors for reduced matrix elements.
 *        Used for operators where tensor order changes between N and overall.
 * @param dir_cos_order Tensor order of the direction cosine part. Renamed ld.
 * @param vib_tensor_order Effective tensor order of the vibrational part or target N-space operator. Renamed lv.
 * @param ixcom_bra Quantum number block for the bra state.
 * @param jxcom_ket Quantum number block for the ket state.
 * @return double The N-dependent correction factor. Returns 0 if selection rules violated.
 */
double rmatrx(const int dir_cos_order, const int vib_tensor_order,
              const int *ixcom_bra, const int *jxcom_ket) /* Renamed ld, lv, ixcom, jxcom */
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
int symnsq(struct SpinvContext *ctx, const int n_sq_power, const int n_s_power,
           const int *iscom_bra_qns, const int *jscom_ket_qns, double *matrix_element_factor)
{
  double n_n_plus_1_power_factor, dot_N_S_bra, dot_N_S_ket, factor_bra_side, factor_ket_side;          /* Renamed x, dotns, dotnsp, zr, zl */
  int delta_2N, first_spin_qn_bra_times_2, N_bra_times_2, N_ket_times_2, J_bra_times_2, J_ket_times_2; /* Renamed modf, ispn, nni, nnj, j */
  unsigned int current_power;                                                                          /* Renamed ipwr */

  N_bra_times_2 = iscom_bra_qns[ctx->nspin];
  N_ket_times_2 = jscom_ket_qns[ctx->nspin];
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
    first_spin_qn_bra_times_2 = iscom_bra_qns[ctx->nspin + 1];                                                                       /* 2*S_bra */
    n_n_plus_1_power_factor = (double)(first_spin_qn_bra_times_2 * (first_spin_qn_bra_times_2 + 2));                                 /* 4*S(S+1) for bra */
    dot_N_S_bra = 0.125 * ((J_bra_times_2 - N_bra_times_2) * (double)(J_bra_times_2 + N_bra_times_2 + 2) - n_n_plus_1_power_factor); /* (J(J+1) - N(N+1) - S(S+1))/2 for bra, times 1/4 from 2J etc. */
                                                                                                                                     /* Actual formula: (J(J+1) - N(N+1) -S(S+1))/2. Here terms are 2X, so ( (2J)(2J+2)/4 - (2N)(2N+2)/4 - (2S)(2S+2)/4 )/2 */
                                                                                                                                     /* = ( J(J+1) - N(N+1) - S(S+1) ) / 2. The 0.125 factor comes from (2X)(2X+2) = 4X(X+1). So (4J(J+1) - 4N(N+1) - 4S(S+1))/8. */

    J_ket_times_2 = jscom_ket_qns[0];                                                                /* 2*J_ket */
    first_spin_qn_bra_times_2 = jscom_ket_qns[ctx->nspin + 1];                                       /* 2*S_ket (reused variable name is for S_ket here) */
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
int symksq(const int k_sq_power, const int K_bra_start, const int K_ket_start,
           const int num_elements, double *work_array_elements,
           short *row_indices, short *col_indices) /* Renamed parameters */
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
int dpmake(const int nsize, double *dp, const double *t, const int n,
           const double *wk, const short *idx, const short *jdx, const int isunit)
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
    memset(dp, 0, nsize * sizeof(double));
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
