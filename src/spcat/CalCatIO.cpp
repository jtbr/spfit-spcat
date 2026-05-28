/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "CalCatIO.hpp"
#include "api/InputSchema.hpp"
#include "splib/calpgm_types.h"
#include "splib/ulib.h"
#include "common/CalError.hpp"
#include "common/file_helpers.hpp"

#define NCARD 130
#define NDVEC 10

bool CalCatIO::readInput(const std::string &intFile,
                         const std::string &varFile,
                         CalCatInput &input,
                         std::unique_ptr<CalculationEngine> &calc_engine,
                         OutputSink *luout)
{
  char titl[NCARD];
  double *dvec;
  int i, j, k, jj, ibcd;
  size_t nl;

  if (!calc_engine) {
    puts("Error: calc_engine is NULL in CalCatIO::readInput.");
    return false;
  }
  if (!luout) {
    puts("Error: luout sink is NULL.");
    return false;
  }

  /* ---- Read .int file ---- */
  FILE *luint = file_helpers::open_input_optional(intFile.c_str());

  /* count lines in .int file */
  for (i = 0; i >= 0; i++) {
    if (fgetstr(titl, NCARD, luint) <= 0)
      break;
  }
  input.ndip = i - 2; /* all but first 2 lines are dipoles */
  rewind(luint);

  /* read first two lines of .int file */
  BOOL first = (fgetstr(titl, NCARD, luint) <= 0);
  if (!first) {
    chtime(titl, 82);
    fputs(titl, stdout);
    luout->puts(titl);
    input.title = std::string(titl);
    first = (fgetstr(titl, NCARD, luint) <= 0);
  }

  dvec = (double *)calalloc((size_t)NDVEC * sizeof(double));
  dvec[0] = 0;
  dvec[1] = 999;
  dvec[2] = 1000;
  dvec[3] = 0;
  dvec[4] = 0;
  dvec[5] = -100;
  dvec[6] = -100;
  dvec[7] = 9999.99;
  dvec[8] = 300.;
  dvec[9] = -1;
  if (!first && pcard(titl, dvec, NDVEC, NULL) == 0)
    first = TRUE;

  input.iflg = (int)dvec[0];
  input.itag = (long)dvec[1];
  input.qrot = dvec[2];
  input.inblk = (int)dvec[3];
  input.lblk = (int)dvec[4];
  input.thrsh = dvec[5];
  input.thrsh1 = dvec[6];
  input.fqmax = dvec[7];
  input.tmq = dvec[8];
  input.maxv = (int)dvec[9];

  if (first) {
    puts(" ERROR IN READING  FIRST CARDS OF .INT");
    fclose(luint);
    free(dvec);
    return false;
  }
  if (input.qrot < 1)
    input.qrot = 1;
  luout->printf("ID=%6ld QSPINROT= %14.4f MIN, MAX QN= %3d %3d\n",
               input.itag, input.qrot, input.inblk, input.lblk);
  luout->printf("MIN.LOG.STR= %9.3f MIN.LOG.STR(FRQ/300GHZ)**2= %9.3f",
               input.thrsh, input.thrsh1);
  luout->printf(" MAX FREQ %10.1f GHZ, TEMP= %8.2f\n", input.fqmax, input.tmq);

  /* read dipole moments */
  int ndip = input.ndip;
  input.dip.resize(ndip);
  input.dip[0] = 0.;
  nl = (size_t)(ndip * sizeof(int));
  input.nvdip.resize(ndip);
  input.nvdip[0] = 0;
  input.isimag.resize(ndip);
  input.isimag[0] = -1;
  nl = (size_t)(ndip * NDECDIP);
  input.idip.resize(nl);
  input.idip[0] = (bcd_t)NDECDIP;

  k = -1;
  ibcd = 0;
  for (j = 0; j < ndip; ++j) {
    input.nvdip[j] = 0;
    input.isimag[j] = -1;
    if (fgetstr(titl, NCARD, luint) <= 0)
      break;
    jj = getbcd(titl, &input.idip[ibcd], NDECDIP);
    if (jj <= 0)
      break;
    dvec[0] = 0.;
    if (pcard(&titl[jj], dvec, 1, NULL) <= 0)
      break;
    input.dip[j] = dvec[0];
    if (NEGBCD(input.idip[ibcd]) == 0 || j == 0) {
      input.nvdip[++k] = 1;
    } else {
      input.nvdip[k] += 1;
    }
    ibcd += NDECDIP;
    luout->printf("%6d %s\n", k + 1, titl);
  }
  input.ndip = j;
  ndip = j;

  input.npdip = 1;
  if (input.iflg >= 10 && (input.iflg % 100) >= 10) {
    /* prstr case: determine npdip from flag */
    if ((input.iflg % 100) >= 20) {
      input.npdip = -(k + 1); /* negative signals "all components" */
    } else {
      input.npdip = 1;
    }
  }
  /* Finalize npdip: if positive, lump all dipoles; if negative, use k+1 */
  if (input.npdip > 0) {
    input.nvdip[0] = ndip;
    input.npdip = 1;
  } else {
    input.npdip = k + 1;
    if (input.npdip > NDVEC) {
      free(dvec);
      dvec = (double *)calalloc((size_t)input.npdip * sizeof(double));
      dvec[0] = 0.;
    }
  }

  fclose(luint); /* close .INT file */

  /* ---- Read .var file ---- */
  FILE *luvar = file_helpers::open_input_optional(varFile.c_str());
  if (!luvar) {
    printf("Error: cannot open .VAR file '%s'\n", varFile.c_str());
    free(dvec);
    return false;
  }

  input.npar = 0;
  input.catqn = MAXCAT;
  if (fgetstr(titl, NCARD, luvar) > 0) {
    luout->puts(".VAR FILE TITLE :");
    luout->puts(titl);
    luout->puts("\n");
    if (fgetstr(titl, NCARD, luvar) > 0) {
      dvec[0] = 0;
      dvec[1] = 1;
      pcard(titl, dvec, 2, NULL);
      input.npar = (int)dvec[0];
      if (dvec[1] < 0.)
        input.catqn = MAXQN;
    }
  }

  /* read option line(s) */
  input.nfmt = input.catqn;
  if (input.npar > 0 && calc_engine->setopt(luvar, &input.nfmt, &input.itd, &input.ndbcd, titl) <= 0)
    input.npar = 0;
  if (input.npar <= 0) {
    puts("Error reading .VAR file");
    if (luvar) fclose(luvar);
    free(dvec);
    return false;
  }

  /* setfmt */
  input.iqnfmtv.resize(input.nfmt << 1, 0);
  input.nqn = calc_engine->setfmt(input.iqnfmtv.data(), input.nfmt);

  /* read parameters and variance */
  int npar = input.npar;
  input.par.resize(npar);
  input.derv.resize(npar);
  nl = (size_t)(npar * input.ndbcd);
  input.idpar.resize(nl);
  input.idpar[0] = (bcd_t)input.ndbcd;
  strcpy(titl, "0123456789");
  int itmp = getpar(luvar, luout->file(), &input.nfit, &npar, input.idpar.data(),
                    input.par.data(), input.derv.data(), titl, 0);
  input.npar = npar;

  unsigned int ndel = (unsigned)(input.nfit + 1);
  size_t nlsq;
  if ((input.nfit & 1) == 0) {
    nlsq = (size_t)((unsigned)input.nfit >> 1) * sizeof(double);
    nlsq *= (size_t)ndel;
  } else {
    nlsq = (size_t)(ndel >> 1) * sizeof(double);
    nlsq *= (size_t)input.nfit;
  }
  input.var.resize(nlsq / sizeof(double));
  getvar(luvar, input.nfit, input.var.data(), input.idpar.data(),
         input.derv.data(), itmp);
  fclose(luvar); /* close .VAR file */

  free(dvec);
  return true;
}

void CalCatIO::write_cat_preamble(FILE *luout, const CatInput &ci, const CalCatInput &cci)
{
    char titl[NCARD];

    /* 1. Title line (from .int/dipoles.toml) + timestamp */
    const std::string &int_t = ci.int_title.empty() ? ci.title : ci.int_title;
    strncpy(titl, int_t.c_str(), NCARD - 27);
    titl[NCARD - 27] = '\0';
    chtime(titl, 82);
    fputs(titl, luout);

    /* 2. ID / QSPINROT / QN line */
    fprintf(luout, "ID=%6ld QSPINROT= %14.4f MIN, MAX QN= %3d %3d\n",
            cci.itag, cci.qrot, cci.inblk, cci.lblk);

    /* 3. LOG.STR / FREQ / TEMP line */
    fprintf(luout, "MIN.LOG.STR= %9.3f MIN.LOG.STR(FRQ/300GHZ)**2= %9.3f",
            cci.thrsh, cci.thrsh1);
    fprintf(luout, " MAX FREQ %10.1f GHZ, TEMP= %8.2f\n", cci.fqmax, cci.tmq);

    /* 4. Dipole lines — reconstruct the .int line format from the BCD id + value */
    {
        int nd = NDECDIP + NDECDIP;
        char sbcd[NDECDIP * 2 + 2];
        int k = -1;
        int ndip = cci.ndip;
        for (int j = 0; j < ndip && j < (int)ci.dipoles.size(); ++j) {
            if (NEGBCD(cci.idip[(size_t)j * NDECDIP]) == 0 || j == 0)
                ++k;
            putbcd(sbcd, nd, &cci.idip[(size_t)j * NDECDIP]);
            const std::string &lbl = ci.dipoles[j].label;
            if (!lbl.empty())
                fprintf(luout, "%6d %s  %.7g / %s\n", k + 1, sbcd, ci.dipoles[j].value, lbl.c_str());
            else
                fprintf(luout, "%6d %s  %.7g\n", k + 1, sbcd, ci.dipoles[j].value);
        }
    }

    /* 5. .VAR FILE TITLE line (from .var/fitted.toml) + timestamp */
    strncpy(titl, ci.title.c_str(), NCARD - 27);
    titl[NCARD - 27] = '\0';
    chtime(titl, 82);
    fputs(".VAR FILE TITLE :", luout);
    fputs(titl, luout);

    /* 6. PARAMETERS - A.PRIORI ERROR header */
    for (int i = 0; i < 30; ++i) fputc(' ', luout);
    fputs("PARAMETERS - A.PRIORI ERROR \n", luout);

    /* 7. One row per parameter */
    {
        int ndbcd = cci.ndbcd;
        int nd = ndbcd + ndbcd;
        char sbcd[32];
        int kfit = 0;
        const double ermin = 1.0e-37;
        double parbase = 1.0;
        int npar = cci.npar;

        for (int i = 0; i < npar; ++i) {
            const bcd_t *id_slot = cci.idpar.data() + (size_t)i * ndbcd;
            bool is_dep = (NEGBCD(id_slot[0]) != 0);
            putbcd(sbcd, nd, id_slot);

            double val = cci.par[i];   /* ratio for dep, absolute for indep */
            double err = cci.derv[i];

            const char *label = "";
            if (i < (int)ci.parameters.size())
                label = ci.parameters[i].label.c_str();

            if (!is_dep) {
                ++kfit;
                parbase = val;
                if (fabs(parbase) < ermin) parbase = ermin;
                fprintf(luout, "%6d %6d %s %21.13E %15.6E %s\n",
                        i + 1, kfit, sbcd, val, err, label);
            } else {
                double abs_val = val * parbase;
                fprintf(luout, "%6d %6d %s %21.13E %15.6f %s\n",
                        i + 1, kfit, sbcd, abs_val, val, label);
            }
        }
        fprintf(luout, "%d parameters read, %d independent parameters\n", npar, kfit);
    }
}
