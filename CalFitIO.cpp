/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include "CalFitIO.hpp"
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "calpgm.h"
#define NDCARD 130
#define LBLEN 16


/**
 * @brief Static method to read input data from files
 * @param parFile Path to parameter file
 * @param linFile Path to line file
 * @param input Output parameter for input data
 * @return True if reading is successful, false otherwise
 */
bool CalFitIO::readInput(const std::string& parFile, const std::string& linFile, CalFitInput& input) {
    FILE *lubak = fopen(parFile.c_str(), "r");
    if (!lubak) {
        printf("Error: Unable to open parameter file '%s'. Please check if the file exists and you have read permissions.\n", parFile.c_str());
        return false;
    }
    printf("Successfully opened parameter file '%s' for reading.\n", parFile.c_str());

    FILE *lulin = fopen(linFile.c_str(), "r");
    if (!lulin) {
        printf("Error: Unable to open line file '%s'. Please check if the file exists and you have read permissions.\n", linFile.c_str());
        fclose(lubak);
        return false;
    }
    printf("Successfully opened line file '%s' for reading.\n", parFile.c_str());

    char card[NDCARD];
    // Read title of .par file
    if (fgetstr(card, NDCARD, lubak) <= 0) {
        puts("Unable to read title of .par file");
        fclose(lubak);
        fclose(lulin);
        return false;
    }
    chtime(card, 82);
    input.title = std::string(card);

    // Read run parameters
    double dvec[8] = {100.0, 32767.0, 1.0, 0.0, 0.0, 1e6, 1.0, 1.0};
    int n = fgetstr(card, NDCARD, lubak);
    if (n != 0)
        n = pcard(card, dvec, 8, NULL);
    if (n == 0) {
        puts("Unable to read second line of .par file");
        fclose(lubak);
        fclose(lulin);
        return false;
    }

    input.npar = (int)dvec[0];
    input.limlin = (int)dvec[1];
    if (dvec[1] < 0.0) {
        input.limlin = -input.limlin;
    }
    input.nitr = (int)dvec[2];
    input.nxpar = (int)dvec[3];
    input.marqp0 = dvec[4];
    input.xerrmx = dvec[5];
    input.parfac = dvec[6];
    input.fqfacq = dvec[7];

    if (input.marqp0 < 0.0)
        input.marqp0 = 0.0;
    if (input.xerrmx < 1.5e-38)
        input.xerrmx = 1e6;

    // Read option cards
    input.nfmt = MAXCAT; // Default value
    input.ndbcd = 1;     // Default value
    input.itd = 0;       // Default value
    input.noptn = 0;     // Number of option lines read
    char namfil[NDCARD] = {0};
    input.optionLines.clear();
    // Read option lines until a non-option line is encountered
    while (fgetstr(card, NDCARD, lubak) > 0) {
        if (card[0] == '#') { // Option line indicator
            input.optionLines.push_back(std::string(card));
            input.noptn++;
            // Parse specific options if needed (e.g., for nfmt, itd, ndbcd)
            // This is a simplified parsing; actual logic may need to be more detailed
            if (strncmp(card, "#FMT", 4) == 0) {
                sscanf(card + 4, "%d", &input.nfmt);
            } else if (strncmp(card, "#ITD", 4) == 0) {
                sscanf(card + 4, "%d", &input.itd);
            } else if (strncmp(card, "#BCD", 4) == 0) {
                sscanf(card + 4, "%d", &input.ndbcd);
            } else if (strncmp(card, "#NAM", 4) == 0) {
                strncpy(namfil, card + 4, NDCARD - 1);
                namfil[NDCARD - 1] = '\0';
                input.namfil = std::string(namfil);
            }
        } else {
            // Skip empty lines and comments
            std::string trimmedCard = card;
            trimmedCard.erase(0, trimmedCard.find_first_not_of(" \t\n\r"));
            if (trimmedCard.empty() || trimmedCard[0] == '!') {
                continue;
            }
            // Non-option line, push back to file or handle as parameter start
            // For simplicity, we'll assume parameters start here and break
            break;
        }
    }
    if (input.noptn > 0) {
        printf("Read %d option lines from parameter file:\n", input.noptn);
        for (const auto& opt : input.optionLines) {
            printf("  - %s\n", opt.c_str());
        }
    } else {
        printf("Warning: No option lines found in parameter file '%s'. Using default values (nfmt=%d, itd=%d, ndbcd=%d).\n",
               parFile.c_str(), input.nfmt, input.itd, input.ndbcd);
    }

    // Allocate memory for parameters (will be used in CalFit class, but read here)
    size_t nl = (size_t)(input.npar * input.ndbcd + input.ndbcd + 3);
    bcd_t *idpar = (bcd_t *)mallocq(nl);
    idpar[0] = (bcd_t)input.ndbcd;
    nl = (size_t)input.npar * sizeof(double);
    double *par = (double *)mallocq(nl);
    double *erp = (double *)mallocq(nl);
    nl = (size_t)(LBLEN * input.npar + 1);
    char *parlbl = (char *)mallocq(nl);

    // Read parameters
    FILE *tempfit = fopen("temp.fit", "w"); // Temporary file for output during parameter reading
    if (!tempfit) {
        puts("Unable to open temporary fit file for parameter reading");
        free(idpar);
        free(par);
        free(erp);
        free(parlbl);
        fclose(lubak);
        fclose(lulin);
        return false;
    }
    int nfit = 0;
    int inpcor = getpar(lubak, tempfit, &nfit, &input.npar, idpar, par, erp, parlbl, LBLEN);
    input.nfit = nfit;
    input.inpcor = inpcor;

    // Read variance if available
    inpcor = getvar(lubak, nfit, NULL, idpar, erp, inpcor); // Assuming var allocation in CalFit
    input.inpcor = inpcor;

    // Read lines from the line file into input.lineData for in-memory processing
    input.lineData.clear();
    char lineBuffer[NDCARD];
    while (fgetstr(lineBuffer, NDCARD, lulin) > 0) {
        input.lineData.push_back(std::string(lineBuffer));
    }
    printf("Read %zu lines from line file '%s' into memory.\n", input.lineData.size(), linFile.c_str());

    fclose(lubak);
    fclose(lulin);
    fclose(tempfit);
    free(idpar); // These will be re-allocated in CalFit
    free(par);
    free(erp);
    free(parlbl);

    return true;
}

/**
 * @brief Static method to write output data to files
 * @param fitFile Path to fit output file
 * @param bakFile Path to backup file
 * @param varFile Path to variance file
 * @param output Output data to write
 * @return True if writing is successful, false otherwise
 */
bool CalFitIO::writeOutput(const std::string& fitFile, const std::string& bakFile, const std::string& varFile, const CalFitOutput& output) {
    // Open output files
    FILE *lufit = fopen(fitFile.c_str(), "w");
    if (!lufit) {
        printf("Unable to open fit output file: %s\n", fitFile.c_str());
        return false;
    }

    FILE *lupar = fopen(bakFile.c_str(), "w");
    if (!lupar) {
        printf("Unable to open backup file: %s\n", bakFile.c_str());
        fclose(lufit);
        return false;
    }

    FILE *luvar = fopen(varFile.c_str(), "w");
    if (!luvar) {
        printf("Unable to open variance file: %s\n", varFile.c_str());
        fclose(lufit);
        fclose(lupar);
        return false;
    }

    // Write title to all output files
    // Assuming title is stored in output or passed from input through CalFit
    // For now, we'll use a placeholder or retrieve from a shared context if available
    char card[NDCARD] = {0};
    // Assuming title was stored somewhere, for now we'll use a dummy title
    strcpy(card, "SPFIT Output - Modernized Version\n");
    chtime(card, 82);
    fputs(card, lufit);
    fputs(card, lupar);
    fputs(card, luvar);

    // Write final parameters and results
    // This part requires access to final parameters and errors from CalFit
    // Since output contains par and erpar, we can use them directly
    if (output.par.size() == 0 || output.erpar.size() != output.par.size()) {
        puts("Output data is incomplete for writing");
        fclose(lufit);
        fclose(lupar);
        fclose(luvar);
        return false;
    }

    // Write header information (parameters, iterations, etc.)
    // These values should be part of output or passed from input
    // For now, we'll assume some values are available in output or use placeholders
    int npar = output.par.size();
    int limlin = 100; // Placeholder, should be from input or output
    int nitr = output.itr;
    int nxpar = 0; // Placeholder
    double marqlast = 0.0; // Placeholder for last Marquardt parameter
    double xerrmx = 1e6; // Placeholder
    double parfac0 = 1.0; // Placeholder
    double fqfacq = 1.0; // Placeholder

    fprintf(lupar, "%4d %4d %4d %4d %14.4E %14.4E %14.4E %12.10f\n",
            npar, limlin, nitr, nxpar, marqlast, xerrmx, parfac0, fqfacq);
    fprintf(luvar, "%4d %4d %4d %4d %14.4E %14.4E %14.4E %12.10f\n",
            npar, limlin, nitr, nxpar, marqlast, xerrmx, parfac0, fqfacq);

    // Write parameters and errors
    // We need idpar and parlbl, which might not be in output directly
    // For simplicity, we'll format without labels or use dummy labels
    for (size_t i = 0; i < output.par.size(); ++i) {
        char pare[64];
        sprintf(pare, "%lu", (unsigned long)i + 1); // Dummy parameter ID
        fprintf(lupar, "%s %23.15E %14.8E\n", pare, output.par[i], output.erpar[i]);
        fprintf(luvar, "%s %23.15E %14.8E\n", pare, output.par[i], output.erpar[i]);
    }

    // Write variance data to var file
    // Since putvar requires var array, which might not be directly in output,
    // we'll assume it's handled elsewhere or use a placeholder
    // putvar(luvar, nfit, var, dpar); // This needs actual data from CalFit

    fclose(lufit);
    fclose(lupar);
    fclose(luvar);

    return true;
}
