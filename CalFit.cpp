/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include <cstring> // for memset
#include <cmath>   // for fabs
#include "lsqfit.h" // for linear algebra functions and other definitions
#include <stdio.h>  // for FILE and fprintf
#include <string.h> // for strcpy
#include <vector>   // for std::vector
#include <algorithm>// for std::copy
#include "CalFitIO.hpp"
#include "CalFit.hpp"
#include "SpinvEngine.hpp"
#include "DpiEngine.hpp"

#define PR_DELAY 6    /* seconds delay between informational messages */
#define FMT_xbgnMW "%5d: %s%14.5f%14.5f %10.5f %10.5f %9.5f"
#define FMT_xblnMW "%14.5f %10.5f %6.4f\n"
#define FMT_xbgnIR "%5d: %s%14.7f%14.7f %10.7f %10.7f %9.7f"
#define FMT_xblnIR "%14.7f %10.7f %6.4f\n"
#define NDCARD 130

/**
 * @brief Constructor for CalFit
 * @param engineType Type of calculation engine to use ("spinv" or "dpi")
 */
CalFit::CalFit(const std::string& engineType) {
    if (engineType == "dpi") {
        calc = std::make_unique<DpiEngine>();
    } else {
        calc = std::make_unique<SpinvEngine>();
    }
}

/**
 * @brief Destructor for CalFit
 */
CalFit::~CalFit() {
    // Unique_ptr will automatically delete the calculation engine
}

/**
 * @brief Run the fitting process
 * @param input Input data for fitting
 * @param output Output data for results
 * @return True if fitting is successful, false otherwise
 */
bool CalFit::run(const CalFitInput& input, CalFitOutput& output) {
    // Initialize parameters using input data from CalFitIO
    if (!initializeParameters(input)) {
        return false;
    }

    // Process spectral lines using input data
    if (!processLines(input)) {
        return false;
    }

    // Perform iterative fitting
    if (!performIteration(input, output)) {
        return false;
    }

    // Calculate final errors and statistics
    if (!calculateErrors(input, output)) {
        return false;
    }

    return true;
}

/**
 * @brief Initialize parameters for fitting
 * @param input Input data for fitting
 * @return True if initialization is successful, false otherwise
 */
bool CalFit::initializeParameters(const CalFitInput& input) {
    // Allocate memory for parameters and related arrays using input data from CalFitIO
    size_t nl = (size_t)(LBLEN * input.npar + 1);
    parlbl = (char *)mallocq(nl);
    nl = (size_t)input.npar * sizeof(double);
    par = (double *)mallocq(nl);
    erp = (double *)mallocq(nl);
    nl = (size_t)(input.npar * input.ndbcd + input.ndbcd + 3);
    idpar = (bcd_t *)mallocq(nl);
    idpar[0] = (bcd_t)input.ndbcd;

    // Set initial values from input
    inpcor = input.inpcor;
    nfit = input.nfit;

    // Placeholder: In a complete implementation, output would be handled via CalFitIO
    // For now, use a temporary file for fit results
    lufit = fopen("temp.fit", "w");
    if (!lufit) {
        puts("Unable to open fit output file");
        return false;
    }

    nl = (size_t)input.npar * sizeof(double);
    oldpar = (double *)mallocq(nl);
    erpar = (double *)mallocq(nl);
    ndfit = nfit + 1;
    ndiag = ndfit + 1;

    // Check if the number of independent parameters fits within memory constraints
    nsize_p = maxmem(&nl); // Assuming maxmem is defined elsewhere or needs to be set
    if (ndfit > nsize_p) {
        printf("number of independent parameters is too big: %d %d\n", nfit, nsize_p - 1);
        return false;
    }

    nl = (size_t)nfit * sizeof(int);
    iperm = (int *)mallocq(nl);
    nl = (size_t)ndfit * sizeof(double);
    dpar = (double *)mallocq(nl);
    delbgn = (double *)mallocq(nl);

    size_t nlsq;
    if ((nfit & 1) == 0) {
        nlsq = (size_t)((unsigned)nfit >> 1) * sizeof(double);
        nlsq *= (size_t)ndfit;
    } else {
        nlsq = (size_t)((unsigned)ndfit >> 1) * sizeof(double);
        nlsq *= (size_t)nfit;
    }
    fitbgn = (double *)mallocq(nlsq);
    var = (double *)mallocq(nlsq);
    nl = (size_t)nfit * sizeof(double);
    nlsq += nl;
    oldfit = (double *)mallocq(nlsq);
    nl *= (size_t)ndfit;
    fit = (double *)mallocq(nl);

    // Fixup supplied values from input
    nxpar = input.nxpar;
    if (nxpar > input.npar || nxpar < 0)
        nxpar = 0;
    if (nxpar > 0) {
        fprintf(lufit, "NUMBER OF PARAMETERS EXCLUDED IF INDEX < 0 =%5d\n", nxpar);
    }
    nxpar = input.npar - nxpar;
    nxfit = nfit;
    for (int i = input.npar - 1; i >= nxpar; --i) {
        int ibcd = i * input.ndbcd;
        if (NEGBCD(idpar[ibcd]) == 0)
            --nxfit;
    }
    xerrmx = input.xerrmx;
    if (xerrmx < tiny)
        xerrmx = 1e6;
    if (fabs(input.fqfacq - 1.) >= 1e-10) {
        fprintf(lufit, " IR Frequencies Scaled by %12.10f\n", input.fqfacq);
    }
    ndfree0 = 0;

    if (inpcor != 0) { // Initialize fit from supplied variance
        double *pvar = var;
        double *pfit = fit;
        for (int n = 1; n <= nfit; ++n) { // Copy from var
            dcopy(n, pvar, 1, pfit, ndfit);
            pvar += n; ++pfit;
        }
        fit[0] = var[0]; dpar[0] = 1.;
        pfit = fit;
        for (int i = 0; i < nfit; ++i) { // Scale
            if (i != 0) {
                pfit += ndfit;
                memset(pfit, 0, sizeof(double) * i);
            }
            int n = nfit - i;
            double *pfitd = pfit + i;
            double val = dnrm2(n, pfitd, 1);
            if (val < tiny) {
                dpar[i] = 1.; scale = 0.;
            } else {
                dpar[i] = scale = 1. / val;
            }
            dscal(n, scale, pfitd, 1);
        }
        pfit = fit;
        for (int i = 0; i < nfit; ++i) {
            int n = i + 1;
            erpar[i] = ddot(n, pfit, ndfit, pfit, ndfit);
            ++pfit;
        }
        // Invert lower triangle
        int n = dqrfac(fit, ndfit, nfit, nfit, 0, erpar, iperm);
        if (nfit > n) {
            fputs("supplied variance is singular\n", lufit);
        }
        dqrsolv(fit, ndfit, n, nfit, 0, iperm);
        pfit = fit;
        for (int n = 1; n <= nfit; ++n) { // Unscale
            scale = dpar[n - 1];
            dscal(n, scale, pfit, ndfit);
            ++pfit;
        }
        fitbgn[0] = fit[0];
        pfit = fit;
        double *pfitb = fitbgn;
        for (int i = 0; i < nfit; ++i) { // Store in fitbgn
            if (i != 0)
                pfit += ndiag;
            int n = nfit - i;
            dcopy(n, pfit, ndfit, pfitb, 1);
            pfitb += n;
        }
    } else { // Initialize fit when there is no supplied variance
        double *pvar = var;
        double *pfitb = fitbgn;
        fitbgn[0] = 1. / var[0];
        var[0] = 0.;
        for (int n = 1; n < nfit; ++n) {
            ++pfitb;
            memset(pfitb, 0, sizeof(double) * n);
            pvar += (n + 1);
            pfitb += n;
            *pfitb = 1. / (*pvar);
            if (*pfitb < 1.e-15)
                --ndfree0;
            *pvar = 0.;
        }
    }
    oldfit[0] = fitbgn[0];
    oldpar[0] = par[0];
    memset(delbgn, 0, sizeof(double) * nfit);
    lbufof(-nfit, input.limlin);

    return true;
}

/**
 * @brief Process spectral lines for fitting
 * @param input Input data for fitting
 * @return True if processing is successful, false otherwise
 */
bool CalFit::processLines(const CalFitInput& input) {
    // Since CalFitIO handles file reading, we process lines directly from input data
    // In a complete implementation, line data would be stored in CalFitInput
    // For now, we'll simulate reading from input data until linein is fully adapted
    int nline = input.limlin;
    int iqnfmt = input.nfmt;

    // Instead of using a temporary file, we will eventually process line data from input.lineData
    // TODO: Implement direct data processing from input.lineData vector or similar structure
    // For now, as a transition, we still use a temporary file but acknowledge the need for change
    FILE *luin = fopen("temp.lin", "r"); // Temporary until linein is updated to use input directly
    if (!luin) {
        puts("Unable to open line input file - ensure temp.lin is created or update to use direct input data");
        return false;
    }

    int mxqn = linein(luin, &nline, iqnfmt);
    fclose(luin);

    if (nline <= 0) {
        puts("No lines read from input file");
        return false;
    }

    // Process lines and set up block structure
    int nblkpf = 1; // Number of blocks per F, to be determined or passed as input
    int flg = -1; // Flag for detailed output, adjust based on input or configuration
    int nbad = lineix(lufit, flg, nline, nblkpf, iqnfmt);

    if (nbad > 0) {
        printf("%d bad lines encountered during processing\n", nbad);
    }

    return true;
}

/**
 * @brief Perform iterative fitting process
 * @param input Input data for fitting
 * @param output Output data for results
 * @return True if iteration is successful, false otherwise
 */
bool CalFit::performIteration(const CalFitInput& input, CalFitOutput& output) {
    static double zero = 0.;
    double *egy, *egyder;
    double *pvar, *pfit, *pfitb, *pfitd;
    double dvec[8], fqfac[4], varv[5], adif, cerr, afrq, rerr, xsqbest, cwid;
    double parfac, ex, xsqir, xsqmw, xerr, xfrq, xsqt, scale, avgir, xwid;
    double avgmw, fqfacq, xerrmx, dif, marqp[3], frq, val, xwt, marqlast;
    double bigd, bigf, bigdIR, bigfIR, bigdMW, bigfMW;
    int ifac, iblk, lblk, iflg, line, icnt, nfir;
    int lstf, nitr, i, k, iblnd, lblnd, nline, initl, noptn, ibase, nf, itd, ibcd, nfmt;
    int limlin, iqnfmt, maxdm, maxf, nrj, nqn, itr, npar, nsize, ndbcd;
    int supblnd, catqn, ndfree, ndfree0;
    int lnext = 0, indx = 0, marqflg = 0, nblkpf = 1; // Initialize variables not defined in class
    short qnum[2*MAXQN];
    char ch, pare[64], aqnum[6*MAXQN+2], card[NDCARD];

    bigdIR = 9.9999999;
    bigdMW = 100. * bigdIR;
    bigfIR = 99999.9999999;
    bigfMW = 100. * bigfIR;
    fqfac[0] = 1;
    fqfac[1] = -1;
    marqp[0] = input.marqp0;
    marqp[1] = -1;
    if (marqp[0] < 0.)
        marqp[0] = 0.;
    limlin = nline = input.limlin;
    nitr = input.nitr;
    npar = input.npar;
    parfac = input.parfac;
    fqfacq = input.fqfacq;
    fqfac[2] = fqfacq / 29979.2458;
    fqfac[3] = -fqfac[2];
    xerrmx = input.xerrmx;
    marqlast = marqp[0];
    ndbcd = input.ndbcd;
    nfmt = input.nfmt;
    itd = input.itd;
    catqn = MAXCAT;
    // Removed reference to input.catqn as it doesn't exist in CalFitInput
    if (input.limlin < 0) {
        catqn = MAXQN;
    }

    // Allocate memory for energy, derivatives, and eigenvector
    maxdm = nsize_p; // This should be set based on earlier processing or input
    size_t nlsq = (size_t) maxdm;
    size_t nl = nlsq * sizeof(double);
    nlsq *= nl;
    double *teig = (double *) mallocq(nlsq);
    nl *= (size_t) (nfit + 2);
    double *pmix = (double *) mallocq(nl);
    pmix[0] = 1.;
    egy = pmix + maxdm;
    egyder = egy + maxdm;

    rqexit(1);
    xsqbest = zero;
    itr = 0;
    dpar[0] = 0.;
    if (nitr < 0)
        nitr = -nitr;
    nsize = 0;

    // Start iteration loop
    do {
        k = 0; ibcd = nxpar * ndbcd;
        for (i = nxpar; i < npar; ++i, ibcd += ndbcd) {
            if (NEGBCD(idpar[ibcd]) == 0)
                dpar[k++] = par[i];
        }
        lstf = lblk = lnext = 0;
        getdbk(&lnext, &iblk, &indx, &initl, &ifac);
        do {
            line = getdbk(&lnext, &iblk, &indx, &initl, &ifac);
            if (iblk != lblk) {       /*  get size of block */
                if (rqexit(0) != 0)
                    break;                /*  check operator interrupt */
                calc->getqn(iblk, 0, 0, qnum, &nsize);
                if (nsize == 0)
                    continue;
                lblk = iblk;
                if (nsize > maxdm) {
                    printf("WARNING .. SIZE OF BLOCK %d2 IS %d AND EXCEEDS DIMENSIONS\n",
                           iblk, nsize);
                    return false;
                }
                k = (iblk - 1) / nblkpf; // Use local variable instead of input.nblkpf
                if (lstf != k && caldelay(PR_DELAY) != 0) {
                    printf("Starting Quantum %3d\n", k);
                    fflush(stdout);
                    lstf = k;
                }
                /*  get energies and derivatives */
                egy = pmix + nsize;
                egyder = egy + nsize;
                calc->hamx(iblk, nsize, npar, idpar, par, egy, teig, egyder, pmix, FALSE);
            }
            /* save energies and derivatives for lines in this block */
            dnuadd(nfit, nxfit, initl, indx, ifac, egy, egyder, nsize, line,
                   dpar, fqfac);
        } while (lnext != 0);       /* repeat until no more energies */
        if (lnext != 0)
            break;
        k = (lblk - 1) / nblkpf; // Use local variable instead of input.nblkpf
        printf("Finished Quantum %3d\n", k);
        fflush(stdout);
        /*   initialize least squares matrix */
        xsqir = xsqmw = avgmw = avgir = 0.;
        pfitb = fitbgn;
        pfit = fit;
        pfitd = fit + nfit;
        dcopy(nfit, &zero, 0, pfitd, ndfit);
        for (int n = 1; n <= nfit; ++n) {
            dcopy(n, pfitb, 1, pfit, ndfit);
            val = delbgn[n - 1];
            daxpy(n, val, pfitb, 1, pfitd, ndfit);
            pfitb += n;
            ++pfit;
        }
        xsqt = ddot(nfit, pfitd, ndfit, pfitd, ndfit);
        nf = nrj = nfir = 0;
        icnt = -1;
        for (i = 0; i < 40; ++i)
            fputc(' ', lufit);
        fputs("EXP.FREQ.  -  CALC.FREQ. -   DIFF.  - EXP.ERR.- ", lufit);
        fputs("EST.ERR.-AVG. CALC.FREQ. -  DIFF. - WT.\n", lufit);
        line = 1; ex = 1.;
        do {                        /*   form least squares matrix, LOOP over lines */
            if (icnt <= 0 && caldelay(PR_DELAY) != 0) {
                printf("Fitting Line %d\n", line);
                fflush(stdout);
                icnt += 50;
            }
            lblnd = line;
            scale = rerr = 0;
            iflg = 0; supblnd = 0; xwid = 1.;
            do {                      /*    UNTIL all elements of blend are found */
                i = frqdat(lblnd, &iblnd, &xfrq, &xwt, &xerr, qnum);
                if ((iblnd & 1) != 0) {
                    supblnd = 1; xwid = xerr;
                    val = 1. / (scale * xwid);
                    val = sqrt(1. - val * val); scale *= val;
                    dscal(nfit + 1, val, dpar, 1);
                } else if (i != 0) {
                    ex = fabs(xwt / xerr);
                    scale += ex;
                    /*  accumulate line contributions */
                    val = dnuget(iflg, nfit, ex, lblnd, dpar);
                    rerr = dpar[nfit];
                    ++iflg;
                }
                ++lblnd;
                --icnt;
            } while (iblnd < 0);
            if (iflg == 0) {
                ++line;
                continue;
            }
            /* calculate errors */
            cerr = calerr(nfit, var, dpar) / scale;
            adif = dpar[nfit] / scale;
            afrq = xfrq - adif;
            if (fabs(rerr) < xerrmx) {
                xsqt += rerr * rerr;
                if (xerr < 0.) {
                    avgir += adif;
                    xsqir += adif * adif;
                    ++nfir;
                } else {
                    avgmw += adif;
                    xsqmw += adif * adif;
                }
                /*  rotate line into FIT matrix */
                jelim(fit, dpar, ndfit, nfit, 1);
                ++nf;
            } else {
                fputs(" ***** NEXT LINE NOT USED IN FIT\n", lufit);
                ++nrj; supblnd = 0;
            }
            if (xerr < 0.) {
                bigf = bigfIR;
                bigd = bigdIR;
            } else {
                bigf = bigfMW;
                bigd = bigdMW;
            }
            if (fabs(afrq) > bigf) {
                afrq = (afrq > 0.) ? bigf : -bigf;
            }
            if (fabs(adif) > bigd) {
                adif = (adif > 0.) ? bigd : -bigd;
            }
            if (iflg == 1) {
                qnfmt2(nqn, qnum, aqnum);
                if (xerr < 0.) {
                    fprintf(lufit, FMT_xbgnIR, line, aqnum, xfrq, afrq, adif, xerr,
                            cerr);
                } else {
                    fprintf(lufit, FMT_xbgnMW, line, aqnum, xfrq, afrq, adif, xerr,
                            cerr);
                }
                fputc('\n', lufit);
                lblnd = line + 1;
            } else {
                lblnd = line; cwid = 0.;
                do {                    /* UNTIL all elements of blend are printed */
                    i = frqdat(lblnd, &iblnd, &xfrq, &xwt, &xerr, qnum);
                    if ((iblnd & 1) != 0) {
                        xfrq = 0.; frq = sqrt(cwid);
                        if (supblnd != 0)
                            xsqt += cwid / (xwid * xwid);
                        qnfmt2(0, qnum, aqnum);
                        if (xerr < 0.) {
                            fprintf(lufit, FMT_xbgnIR, lblnd, aqnum, xfrq, frq, frq, xwid,
                                    xfrq);
                        } else {
                            fprintf(lufit, FMT_xbgnMW, lblnd, aqnum, xfrq, frq, frq, xwid,
                                    xfrq);
                        }
                        fputc('\n', lufit);
                    } else if (i != 0) {
                        if (supblnd != 0) {
                            ex = sqrt(xwt) / xwid;
                            frq = dnuget(0, nfit, ex, lblnd, dpar);
                            jelim(fit, dpar, ndfit, nfit, 1);
                            ++nf;
                        } else {
                            frq = dnuget(-1, nfit, ex, lblnd, dpar);
                        }
                        dif = xfrq - frq;
                        if (fabs(frq) > bigf)
                            frq = (frq > 0.) ? bigf : -bigf;
                        if (fabs(dif) > bigd)
                            dif = (dif > 0.) ? bigd : -bigd;
                        qnfmt2(nqn, qnum, aqnum);
                        if (xerr < 0.) {
                            fprintf(lufit, FMT_xbgnIR, lblnd, aqnum, xfrq, frq, dif, xerr,
                                    cerr);
                            fprintf(lufit, FMT_xblnIR, afrq, adif, xwt);
                        } else {
                            fprintf(lufit, FMT_xbgnMW, lblnd, aqnum, xfrq, frq, dif, xerr,
                                    cerr);
                            fprintf(lufit, FMT_xblnMW, afrq, adif, xwt);
                        }
                        dif = frq - afrq;
                        cwid += xwt * dif * dif;
                    }
                    ++lblnd;
                } while (iblnd < 0);
            }
            line = lblnd;
        } while (line <= nline);    /* end loop over lines */
        if (nrj > 0) {
            printf(       "%5d Lines rejected from fit\n", nrj);
            fprintf(lufit,"%5d Lines rejected from fit\n", nrj);
        }
        if (nf < 1)
            nf = 1;
        /* zero upper triangle of fit matrix */
        pfit = fit;
        for (k = 1; k < nfit; ++k) {
            pfit += ndfit;
            memset(pfit, 0, sizeof(double) * k);
        }
        varv[0] = xsqt + nrj * xerrmx * xerrmx;
        marqlast = marqp[0];
        marqflg = lsqfit(fit, ndfit, nfit, 1, marqp, varv,
                         oldfit, erpar, dpar, iperm);
        if (marqflg != 0) {
            dcopy(input.npar, oldpar, 1, par, 1);
            strcpy(card, "Fit Diverging: restore parameters\n");
            fputs(card, lufit);
            fputs(card, stdout);
        } else {
            xsqbest = (xsqt - varv[3]) / nf;
            if (xsqbest > 0.)
                xsqbest = sqrt(xsqbest);
            /* print normalized diagonals */
            fputs("NORMALIZED DIAGONAL:\n", lufit);
            icnt = 0;
            for (k = 0; k < nfit; ++k) {
                fprintf(lufit, "%5d %13.5E", k + 1, erpar[k]);
                if ((++icnt) == 6) {
                    icnt = 0;
                    fputc('\n', lufit);
                }
            }
            if (icnt > 0)
                fputc('\n', lufit);
        }
        printf(       "MARQUARDT PARAMETER = %g, TRUST EXPANSION = %4.2f\n",
                  marqp[0], marqp[2]);
        fprintf(lufit,"MARQUARDT PARAMETER = %g, TRUST EXPANSION = %4.2f\n",
                  marqp[0], marqp[2]);
        if (input.parfac < 0. && xsqt >= 0.) {
            ndfree = nf + ndfree0;
            if (ndfree <= 0)
                ndfree = 1;
            parfac = -input.parfac * sqrt(xsqt / ndfree);
            fputs("WARNING: parameter errors multiplied by ", lufit);
            fprintf(lufit, "%15.5f for %8d degrees of freedom\n",
                    parfac, ndfree);
        }
        xsqt = xsqt / nf;
        if (xsqt > 0.)
            xsqt = sqrt(xsqt);
        if (xsqbest > xsqt)
            xsqbest = xsqt;
        /*   get estimated errors  and print parameters */
        for (i = 0; i < 32; ++i)
            fputc(' ', lufit);
        fputs("NEW PARAMETER (EST. ERROR) -- CHANGE THIS ITERATION\n", lufit);
        pfitd = fit;
        char *tlblnxt = parlbl;
        k = ibase = 0;
        for (i = 0, ibcd = 0; i < npar; ++i, ibcd += ndbcd) {
            char *tlbl = tlblnxt;
            tlblnxt += LBLEN;
            oldpar[i] = par[i];
            if (NEGBCD(idpar[ibcd]) == 0) {
                if (k != 0)
                    pfitd += ndiag;
                int n = nfit - k;
                dif = pfitd[n];
                par[i] += dif;
                delbgn[k] -= dif;
                erpar[i] = dpar[k] = parfac * dnrm2(n, pfitd, 1);
                dscal(n, parfac, pfitd, 1);
                parer(par[i], erpar[i], dif, pare);
                putbcd(card, NDCARD, &idpar[ibcd]);
                ch = (*tlblnxt); *tlblnxt = '\0';
                fprintf(lufit, "%4d %s %10.10s %s\n", ++k, card, tlbl,
                        pare);
                *tlblnxt = ch;
                ibase = i;
            }
        }
        pfit = fit;
        pvar = var;
        for (int n = 1; n <= nfit; ++n) {       /* copy to var */
            dcopy(n, pfit, ndfit, pvar, 1);
            ++pfit;
            pvar += n;
        }
        if (nf > nfir) {
            scale = 1. / (double) (nf - nfir);
            avgmw *= scale;
            xsqmw = sqrt(xsqmw * scale);
        }
        if (nfir > 0) {
            scale = 1. / (double) nfir;
            avgir *= scale;
            xsqir = sqrt(xsqir * scale);
        }
        printf(" MICROWAVE AVG = %15.6f MHz, IR AVG =%15.5f\n",
                      avgmw, avgir);
        fprintf(lufit," MICROWAVE AVG = %15.6f MHz, IR AVG =%15.5f\n",
                      avgmw, avgir);
        printf(       " MICROWAVE RMS = %15.6f MHz, IR RMS =%15.5f\n",
                      xsqmw, xsqir);
        fprintf(lufit," MICROWAVE RMS = %15.6f MHz, IR RMS =%15.5f\n",
                      xsqmw, xsqir);
        ++itr;
        printf(       " END OF ITERATION %2d OLD, NEW RMS ERROR=%15.5f %15.5f\n",
                      itr, xsqt, xsqbest);
        fprintf(lufit," END OF ITERATION %2d OLD, NEW RMS ERROR=%15.5f %15.5f\n",
                      itr, xsqt, xsqbest);
        fflush(stdout);
    } while (itr < nitr && 0.999999 * xsqt > xsqbest);

    // Store results in output
    output.xsqbest = xsqbest;
    output.itr = itr;

    free(pmix);
    free(teig);

    return true;
}

/**
 * @brief Calculate errors and final statistics
 * @param input Input data for fitting to access necessary parameters
 * @param output Output data for results
 * @return True if calculation is successful, false otherwise
 */
bool CalFit::calculateErrors(const CalFitInput& input, CalFitOutput& output) {
    lbufof(-1, 0);      /* release storage */
    int nblkpf = 1, maxdm = nsize_p; // Define variables for setblk call
    calc->setblk(lufit, 0, idpar, par, &nblkpf, &maxdm); /* release storage */
    if (output.itr == 0) {
        puts(" output files not updated");
        return false;
    }
    dcopy(input.npar, oldpar, 1, par, 1);

    /*  compute correlation matrix */
    prcorr(lufit, nfit, fit, ndfit, dpar);
    char card[NDCARD] = {0}; // Define card for output
    fputs(card, lufit);
    fclose(lufit);

    // Update output with final results
    int localNpar = input.npar; // Use input.npar since it's available from input
    output.par.resize(localNpar);
    output.erpar.resize(localNpar);
    std::copy(par, par + localNpar, output.par.begin());
    std::copy(erpar, erpar + localNpar, output.erpar.begin());

    free(fit);
    free(oldfit);
    free(var);
    free(fitbgn);
    free(delbgn);
    free(dpar);
    free(iperm);
    free(erpar);
    free(oldpar);
    free(erp);
    free(idpar);
    free(par);
    free(parlbl);

    return true;
}

/**
 * @brief Format quantum numbers as a string for output
 * @param nqn Number of quantum numbers to format
 * @param qnum Array of quantum numbers
 * @param aqnum Output string buffer for formatted quantum numbers
 * @return Always returns 0
 */
int CalFit::qnfmt2(int nqn, short *qnum, char *aqnum) {
    // Implementation copied from calfit.cpp
    int i;
    for (i = 0; i < nqn; ++i) {
        sprintf(aqnum, "%3d", (int) qnum[i]);
        aqnum += 3;
    }
    for (i = nqn; i < 12; ++i) {
        aqnum[2] = aqnum[1] = aqnum[0] = ' ';
        aqnum += 3;
    }
    aqnum[0] = '\0';
    return 0;
}

/**
 * @brief Format parameter values, errors, and changes for output
 * @param par Parameter value
 * @param errx Parameter error
 * @param dif Parameter change
 * @param ptmp Output string buffer for formatted parameter
 * @return Always returns 0
 */
int CalFit::parer(double par, double errx, double dif, char *ptmp) {
    static int czero = (int) '0';  /* ASCII code for '0' character */
    char *pfmt;                    /* Pointer for building format string */
    double adif, apar, aten, aerr; /* Working copies of input values */
    char chexp[6], fmt[34];        /* Format string components */
    int msd, id, ie, efield, ip, lsd, k; /* Format calculation variables */

    /*      sets up special format for parameters and errors */
    /*     PAR  = parameter value */
    /*     ERRX = parameter error */
    /*     DIF  = parameter change */
    /*     PTMP = output string for printing */

    apar = par;
    aerr = errx;
    adif = dif;
    efield = 0;
    aten = 1.;
    /*     compute exponent fields */
    ie = (int) (log10(fabs(aerr) + 1.e-37) - 102.5) + 100;
    id = (int) (log10(fabs(adif) + 1.e-37) - 100.0) + 100;
    ip = (int) (log10(fabs(apar) + 1.e-37) - 100.0) + 100;
    lsd = -ie;
    if (lsd < 0)
        lsd = 0;
    msd = (ip > id) ? ip : id;
    /*  check for too many digits */
    k = 14 - ip;
    if (k < lsd)
        lsd = k;
    k = 10 - id;
    if (k < lsd)
        lsd = k;
    if (msd <= -2) {              /* number too small without exponent */
        k = (1 - msd) / 3;
        efield = -3 * k;
        while ((--k) >= 0)
            aten *= 1000;
    } else if (lsd < 0) {         /* number too big without exponent */
        k = (1 + msd) / 3;
        if (k > 0)
            efield = 3 * k;
        while ((--k) >= 0)
            aten *= 0.001;
    }
    if (efield != 0) {            /* E format */
        lsd += efield;
        memcpy(chexp, "0fE+00", 6);
        if (efield < 0) {
            chexp[3] = '-';
            efield = -efield;
        }
        msd = efield / 10;
        if (msd > 0) {
            efield -= msd * 10;
            chexp[4] = (char) (msd + czero);
        }
        chexp[5] = (char) (efield + czero);
        apar *= aten;
        aerr *= aten;
        adif *= aten;
    } else {                      /* F format */
        memcpy(chexp, "0f    ", 6);
    }
    if (lsd > 9)
        lsd = 9;
    if (lsd > 0)
        chexp[0] = (char) (lsd + czero);
    while ((lsd--) > 0)
        aerr *= 10.;
    ie = (int) (aerr + 0.5);
    pfmt = fmt;
    memcpy(pfmt, "%16.", 4);
    pfmt += 4;
    memcpy(pfmt, chexp, 2);
    pfmt += 2;
    memcpy(pfmt, "(%3d)", 5);
    pfmt += 5;
    memcpy(pfmt, chexp + 2, 4);
    pfmt += 4;
    memcpy(pfmt, " %12.", 5);
    pfmt += 5;
    memcpy(pfmt, chexp, 6);
    pfmt += 6;
    *pfmt = '\0';
    sprintf(ptmp, fmt, apar, ie, adif);
    return 0;
}

/**
 * @brief Read experimental spectral lines from input file
 * @param luin Input file pointer
 * @param nline Pointer to number of lines (input: max lines, output: actual lines)
 * @param iqnfmt Quantum number format for line input
 * @return Largest quantum number encountered
 */
int CalFit::linein(FILE *luin, int *nline, int iqnfmt) {
    /* Local variables */
    SXLINE *xline;                /* Pointer to line data structure */
    double xfrqn, xerrn, xwtn;    /* Current line frequency, error, weight */
    double xfrqx, xerrx;          /* Previous line frequency and error */
    int nqn, nqnu, nqnl;          /* Number of quantum numbers */
    int kqnu, kqnl;               /* Indices for upper/lower state quantum numbers */
    int i, iqf, ipace, mxline;    /* Loop variables and counters */
    int mxqn, isblnd, icmp;       /* Max quantum number, blend flag, comparison flag */
    short nbln, nqnt[20], *iqnum; /* Blend counter, quantum number template, quantum number pointer */
    static char card[NDCARD];     /* Buffer for reading input lines */

    /*   get lines from input  and stores them */
    /*     LUIN= unit for finding lines */
    /*     NLINE = number of lines */
    /*     IQNFMT= qunatum number format for line input */
    /*     RETURN: largest quantum number */
    /*******************************************************************/

    mxline = *nline;
    mxqn = 1;
    nbln = 1;
    nqn = deflin(iqnfmt, nqnt);
    nqnu = nqn - 1;
    if (nqnt[nqnu] < 0)
        nqnu = 0;
    kqnu = nqnt[nqnu];
    nqnl = nqnu + nqn;
    kqnl = nqnt[nqnl];
    ipace = 100;
    xfrqx = xerrx = 0.;
    icmp = 0;
    for (i = 1; i <= mxline; ++i) {       /*  loop for reading lines */
        xline = lbufof(1, i);
        iqnum = xline->qn;
        if (getlin(luin, nqn, nqnt, iqnum, &xfrqn, &xerrn, &xwtn,
                   card, NDCARD) < 0) {
            *nline = i - 1;
            return mxqn;
        }
        iqf = iqnum[nqnu];
        if (iqf == -1) {
            if (kqnu >= 0) {
                iqf = -iqnum[kqnu];
                if (iqf >= 0)
                    iqf = -1;
            }
            iqnum[nqnu] = (short) iqf;
        }
        if (iqf < 0)
            iqf = -iqf;
        if (iqf > mxqn)
            mxqn = iqf;
        iqf = iqnum[nqnl];
        if (iqf == -1) {
            if (kqnl > 0) {
                iqf = -iqnum[kqnl];
                if (iqf >= 0)
                    iqf = -1;
            }
            iqnum[nqnl] = (short) iqf;
        }
        if (iqf < 0)
            iqf = -iqf;
        if (mxqn < iqf)
            mxqn = iqf;
        xline->xfrq = xfrqn;
        xline->xerr = (float) xerrn;
        xline->xwt = (float) fabs(xwtn);
        xline->linku = 0;
        xline->linkl = 0;
        isblnd = 0;
        if (icmp != 0 && fabs(xfrqn - xfrqx) < fabs(xfrqn) * 1.e-14 + 1.e-8) {
            /* frq match */
            if (fabs(xerrn - xerrx) < 1e-7) {
                isblnd = 1;
            } else if ((xerrn / xerrx) > 2.0 && nbln > 2) {
                isblnd = 1; ++nbln; icmp = 0;
                xline->xwt = (float)0.;
                iqnum[0] = (short)-1;
                iqnum[nqn] = iqnum[0];
            }
        }
        if (isblnd != 0) {
            xline->bln = nbln;
            xline = lbufof(1, i - 1);
            xline->bln = -2;
            nbln += 2;
        } else {
            xline->bln = 0;
            nbln = 2; icmp = 1;
        }
        if (ipace <= i) {
            ipace += 100;
            printf("Reading Line %d\n", i);
            fflush(stdout);
        }
        xerrx = xerrn;
        xfrqx = xfrqn;
    }
    return mxqn;
}

/**
 * @brief Process spectral lines and set up block structure for fitting
 * @param lu Output file pointer for listing
 * @param flg Flag for printing (negative for detailed output)
 * @param nline Number of lines
 * @param nblkpf Number of blocks per F quantum number
 * @param iqnfmt Quantum number format for line input
 * @return Number of bad lines
 */
int CalFit::lineix(FILE *lu, int flg, int nline, int nblkpf, int iqnfmt) {
    /*   get lines from input and store them */
    /*     LU = unit for printout of lines ( if > 0 ) */
    /*     NLINE = number of lines */
    /*     NBLKPF= number of blocks per F */
    /*     IQNFMT= qunatum number format for line input */
    /******************************************************************/
    static int nsort = 2048;
    SXLINE *xline;
    double xfrqn, xerrn, xwtn, xnorm;
    int nblk, ipos, i, j, ipace, nread, iblkl, iblku, ncat;
    int linkx, indxl, linky, indxu, orgblk, nqn, nqn2, nbad;
    /*@owned@*/ int *prvblk;
    short *iqnum;
    char aqnum[6*MAXQN+2];

    nbad = 0;
    nblk = 0;
    prvblk = (int *) mallocq((size_t) (nsort + 1) * sizeof(int));
    prvblk[0] = 0;
    for (i = 1; i <= nsort; ++i) {
        prvblk[i] = 0;
    }
    nqn = iqnfmt % 10;
    if (nqn == 0) nqn = 10;
    nqn2 = nqn + nqn; ncat = nqn2;
    if (ncat < 12) ncat = 12;
    i = (iqnfmt / 100) % 5;
    if (i >= nqn) {
        ipos = 1;
    } else {
        ipos = nqn;
    }
    if (flg < 0) {
        fputs(" LINE,BLKU,INDXU,BLKL,INDXL,QUANTUM NUMBERS", lu);
        for (i = 0; i < 19; ++i)
            fputc(' ', lu);
        fputs("ENERGY    EXP. ERROR    WEIGHTS\n", lu);
    }
    xnorm = 0.;
    ipace = 50;
    /*       loop for converting lines */
    for (nread = 1; nread <= nline; ++nread) {
        xline = lbufof(1, nread);
        xfrqn = xline->xfrq;
        xerrn = xline->xerr;
        xwtn = xline->xwt;
        /* find blocks and index for upper and lower states */
        iqnum = xline->qn;
        getblk(&iblku, &indxu, iqnum, nblkpf, ipos, nqn);
        getblk(&iblkl, &indxl, &iqnum[nqn], nblkpf, ipos, nqn);
        if (iblkl == 0 && iqnum[nqn] >= 0)
            iblku = 0;
        xline->ibu = iblku;
        xline->inu = (short) indxu;
        xline->ibl = iblkl;
        xline->inl = (short) indxl;
        if (iblku == 0 && (xline->bln & 1) == 0) {
            /*  print out bad line and try for next */
            ++nbad;
            xline->xwt = 0.;
            xwtn = 0.;
            qnfmt2(nqn2, iqnum, aqnum);
            printf(    "Bad Line(%3d): %s %14.5f %8.5f\n",
                       nread, aqnum, xfrqn, xerrn);
            fprintf(lu,"Bad Line(%3d): %s %14.5f %8.5f\n",
                       nread, aqnum, xfrqn, xerrn);
        } else {
            /*  set up links for calculating in order of block */
            if (iblku <= iblkl) {
                lnlink(prvblk, nsort, iblku, nread);
                lnlink(prvblk, nsort, iblkl, -nread);
                if (nblk < iblkl)
                    nblk = iblkl;
            } else {
                lnlink(prvblk, nsort, iblkl, -nread);
                lnlink(prvblk, nsort, iblku, nread);
                if (nblk < iblku)
                    nblk = iblku;
            }
            if (flg < 0) {
                iqnum = xline->qn;
                fprintf(lu," %4d%4d%4d%4d%4d:", nread, iblku, indxu, iblkl, indxl);
                for (i = 0; i < ncat; ++i) {
                    j = iqnum[i];
                    fprintf(lu, "%3d", j);
                }
                fprintf(lu, " %14.4f %9.4f %9.4f", xfrqn, xerrn, xwtn);
                j = xline->bln;
                if (j != 0) {
                    fprintf(lu, "   Line Blended with %3d\n", nread - (j >> 1));
                } else {
                    fputc('\n', lu);
                }
            }
        }
        /* let the user know something is happening */
        if (ipace <= nread || nread == nline) {
            ipace += 50;
            printf("Converting Line %d\n", nread);
            fflush(stdout);
        }
        j = xline->bln;
        if (j != 0) {
            xnorm += xwtn;
            if (j > 0) {     /* normalize weights */
                xnorm = 1. / xnorm;
                for (j = nread - (j >> 1); j <= nread; ++j) {
                    xline = lbufof(1, j);
                    xline->xwt = (float) (xline->xwt * xnorm);
                }
                xnorm = 0.;
            }
        } else {
            xline->xwt = 1.;
        }
    }                             /* end loop for converting lines */
    orgblk = nsort;
    while (nblk > orgblk) {       /* finish up links */
        --orgblk;
        linkx = 0;
        for (j = 0; j < nsort; ++j) {
            if (prvblk[j] != 0) {
                linkx = prvblk[j];
                prvblk[j] = 0;
            }
        }
        prvblk[0] = linkx;
        prvblk[nsort] = 0;
        linkx = lnlink(prvblk, nsort, 1, 0);
        while (linkx != 0) {
            linky = linkx;
            getdbk(&linkx, &iblkl, &j, &j, &j);
            iblkl -= orgblk;
            lnlink(prvblk, nsort, iblkl, linky);
        }
        orgblk += nsort;
    }
    free(prvblk);
    return nbad;
}

/**
 * @brief Get block and index from quantum numbers
 * @param iblk Output block number
 * @param indx Output index within block
 * @param iqnum Array of quantum numbers
 * @param nblkpf Number of blocks per F quantum number
 * @param ipos Position in iqnum to find F
 * @param nqn Number of quantum numbers
 * @return Always returns 0
 */
int CalFit::getblk(int *iblk, int *indx, short *iqnum, int nblkpf, int ipos, int nqn) {
    int ibgn, kbgn, nblk, iqnf, k, idblk, kdtau, kk, nn, icmp;
    int iblkx, indxx, idgn, nsize, ncod;
    short kqnum[MAXQN];

    /*  gets block (IBLK) and INDEX from quantum numbers (IQNUM) */
    /*     NBLKPF IS THE NUMBER OF BLOCKS PER "F" */
    /*     IPOS   IS THE POSITION IN IQNUM TO FIND "F" */
    /*     NQN    IS THE NUMBER OF QUANTUM NUMBERS */

    *indx = 0;
    *iblk = 0;
    /*   ..check for indication of null level */
    if (iqnum[0] < 0)
        return 0;
    --ipos;
    iqnf = iqnum[ipos];
    iqnum[ipos] = iqnum[0];
    ibgn = iqnf;
    nblk = nblkpf;
    if (ibgn < 0)
    {
        ibgn = -ibgn;
        k = 10 - ibgn;
        if (k > 1)
            nblk *= k;
    }
    ibgn *= nblkpf;
    /*   ..loop over blocks for given F */
    icmp = 1;
    iblkx = indxx = 0;
    for (idblk = 1; idblk <= nblk; ++idblk)
    {
        iblkx = ibgn + idblk;
        calc->getqn(iblkx, 0, 0, kqnum, &nsize);
        indxx = 1;
        /*       ..search for match within block */
        while (indxx <= nsize)
        {
            ncod = calc->getqn(iblkx, indxx, nqn, kqnum, &idgn);
            kqnum[ipos] = kqnum[0];
            nn = (ncod > 0) ? ncod : -ncod;
            /*           ..NN is the size of a wang sub-block */
            /*           ..NCOD < 0 if oblate basis */
            if (nn <= 1)
            {
                nn = 1;
                kbgn = 1;
            }
            else
            {
                kbgn = 3;
            }
            icmp = 0;
            for (k = kbgn; k < nqn; ++k)
            {
                icmp = iqnum[k] - kqnum[k];
                if (icmp != 0)
                    break;
            }
            if (icmp == 0)
            {
                /* all quanta tested equal */
                if (nn == 1)
                    break;
                /* assignment could be in this sub-block */
                kdtau = iqnum[1] - iqnum[2] - kqnum[1] + kqnum[2];
                if (ncod < 0)
                    kdtau = -kdtau;
                if (kdtau >= 0 && (kdtau & 3) == 0)
                { /* symmetry is good */
                    kk = (int)((unsigned)kdtau >> 2);
                    if (kk < nn)
                    {
                        /* value is in range */
                        indxx += kk;
                        break;
                    }
                }
                ++icmp;
            }
            /*  to get here, quanta were not right */
            indxx += nn;
        }
        if (icmp == 0)
            break;
    }
    if (icmp == 0)
    { /*   standard return */
        if (nblk > nblkpf)
        {
            iqnf -= (iblkx - ibgn - 1) / nblkpf;
        }
        if (iqnf < 0)
            indxx = -indxx;
        *iblk = iblkx;
        *indx = indxx;
    }
    iqnum[ipos] = (short)iqnf;
    return 0;
}
