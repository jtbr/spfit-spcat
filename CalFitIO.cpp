/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

/*   Herbert M. Pickett, 20 March 1989 */
/*   Revised version in C++, 2025 for modernization */

#include <stdio.h>
#include <string.h>
#include <algorithm> // For std::copy
#include <vector>    // For std::vector
#include "CalFitIO.hpp"
#include "calpgm.h"  // For bcd_t, fgetstr, chtime, pcard, getpar, getvar, mallocq, MAXCAT etc.


/**
 * @brief Static method to read input data from files
 * @param parFile Path to parameter file
 * @param linFile Path to line file
 * @param input Output parameter for input data
 * @param lufit_for_getpar File stream for getpar to write its log.
 * @return True if reading is successful, false otherwise
 */
bool CalFitIO::readInput(const std::string &parFile, const std::string &linFile,
                         CalFitInput &input, FILE *lufit_for_getpar)
{
  if (!lufit_for_getpar)
  {
    puts("Error: lufit_for_getpar stream is NULL in CalFitIO::readInput.");
    return false;
  }

  FILE *lubak = fopen(parFile.c_str(), "r");
  if (!lubak)
  {
    printf("Error: Unable to open parameter file '%s'.\n", parFile.c_str());
    return false;
  }
  printf("Successfully opened parameter file '%s' for reading.\n", parFile.c_str());

  char card[NDCARD];

  // Read title of .par file
  if (fgetstr(card, NDCARD, lubak) <= 0)
  {
    puts("Unable to read title of .par file");
    fclose(lubak);
    return false;
  }
  chtime(card, 82);
  input.title = std::string(card);
  fputs(input.title.c_str(), lufit_for_getpar); // Write title to lufit

  // Read run parameters
  double dvec[8] = {100.0, 32767.0, 1.0, 0.0, 0.0, 1e6, 1.0, 1.0}; // Defaults from original
  int n_read_dvec = fgetstr(card, NDCARD, lubak);
  if (n_read_dvec != 0)
    n_read_dvec = pcard(card, dvec, 8, NULL);
  if (n_read_dvec == 0)
  {
    puts("Unable to read second line of .par file");
    fclose(lubak);
    return false;
  }

  input.npar = (int)dvec[0];
  input.limlin = (int)dvec[1]; // Store original value, sign matters for catqn check later
  // The actual number of lines might change after linein, limlin is the request/limit
  input.nitr = (int)dvec[2];
  input.nxpar_from_file = (int)dvec[3];
  input.marqp0 = dvec[4];
  input.xerrmx = dvec[5];
  input.parfac_initial = dvec[6];
  input.fqfacq = dvec[7];

  // Defaults for options
  input.nfmt_from_options = MAXCAT; // MAXCAT from calpgm.h
  input.ndbcd_from_options = 1;     // Default BCD length if not specified
  input.itd_from_options = 0;       // Default ITD
  input.noptn_count = 0;
  input.namfil_from_options.clear();
  input.optionLines.clear();

  // Read option cards
  // This part is tricky as calc->setopt in original might have complex logic
  // We replicate the simplified parsing for now and store raw lines
  // In CalFit, calc->setopt would be called with these raw lines or a temp file
  // For now, CalFitIO just parses some known ones.
  std::string temp_namfil_str; // Temporary for namfil from options

  // The original setopt call happens in main *before* getpar
  // It's possible setopt influences ndbcd, which is needed for getpar
  // For now, we parse options here to get ndbcd.
  // A more robust solution would be to pass lubak to CalFit for setopt,
  // then CalFit passes back ndbcd, or CalFit does getpar internally.
  // Sticking to current plan: CalFitIO reads options, CalFit uses parsed values.

  long current_pos = ftell(lubak); // Save position before reading options

  while (fgetstr(card, NDCARD, lubak) > 0)
  {
    std::string trimmedCard = card;
    trimmedCard.erase(0, trimmedCard.find_first_not_of(" \t\n\r"));
    trimmedCard.erase(trimmedCard.find_last_not_of(" \t\n\r") + 1);

    if (trimmedCard.empty() || trimmedCard[0] == '!')
    {                             // Skip empty lines and comments
      current_pos = ftell(lubak); // Update position
      continue;
    }

    // Check if it's an option line (e.g., starts with '#' or specific keywords if not '#')
    // The original `setopt` in CalculationEngine would know how to identify option lines.
    // Assuming option lines are marked or `getpar` knows when parameters start.
    // For simplicity, let's assume `getpar` robustly finds the start of parameters.
    // The original `setopt` reads until it finds a non-option line.
    // We need to replicate that behavior. Let's assume options start with '#' for this example.
    // If your options are not prefixed, this logic needs to be more like original setopt.

    bool is_option = false; // This needs to be determined by how setopt works
    // Simplified check:
    if (card[0] == '#')
      is_option = true; // Example: options start with #
    // Or if using SpinvEngine/DpiEngine's setopt logic:
    // You might need a way to "peek" if the line is an option without consuming it
    // from lubak if the line is to be passed to getpar.

    // This is a placeholder for more robust option detection.
    // For now, we assume CalFitIO can parse these specific options if present
    // AND that these option lines won't be mistaken for parameter lines by getpar.
    // It's better if `setopt` itself is used.
    // For Step 2, we will integrate `calc->setopt`.
    // For now, let's parse known options if they appear.

    if (strncmp(card, "#FMT", 4) == 0)
    { // Example simplified parsing
      sscanf(card + 4, "%d", &input.nfmt_from_options);
      is_option = true;
    }
    else if (strncmp(card, "#ITD", 4) == 0)
    {
      sscanf(card + 4, "%d", &input.itd_from_options);
      is_option = true;
    }
    else if (strncmp(card, "#BCD", 4) == 0)
    {
      sscanf(card + 4, "%d", &input.ndbcd_from_options);
      is_option = true;
    }
    else if (strncmp(card, "#NAM", 4) == 0)
    {
      char temp_nam_buf[NDCARD] = {0};
      strncpy(temp_nam_buf, card + 5, NDCARD - 6); // Skip "#NAM "
      temp_nam_buf[NDCARD - 6] = '\0';             // Ensure null termination
      input.namfil_from_options = std::string(temp_nam_buf);
      // Remove trailing newline if present from fgetstr
      input.namfil_from_options.erase(input.namfil_from_options.find_last_not_of(" \n\r\t") + 1);
      is_option = true;
    }
    // Add more specific option parsers as needed from setopt logic

    if (is_option)
    {
      input.optionLines.push_back(std::string(card));
      input.noptn_count++;
      fputs(card, lufit_for_getpar); // Write option line to lufit
      if (card[strlen(card) - 1] != '\n')
        fputc('\n', lufit_for_getpar);
      current_pos = ftell(lubak); // Update position after consuming an option line
    }
    else
    {
      // Not a recognized option line by this simplified parser, assume parameters start
      fseek(lubak, current_pos, SEEK_SET); // Rewind to before this line for getpar
      break;
    }
  }
  // End of simplified option reading

  // Allocate temporary C-style arrays for getpar
  // Size for idpar depends on npar and ndbcd_from_options
  size_t idpar_elem_count = (size_t)input.npar * input.ndbcd_from_options + input.ndbcd_from_options + 3;
  unsigned char *temp_idpar = (unsigned char *)mallocq(idpar_elem_count * sizeof(unsigned char));
  if (!temp_idpar)
  {
    puts("mallocq failed for temp_idpar");
    fclose(lubak);
    return false;
  }
  temp_idpar[0] = (unsigned char)input.ndbcd_from_options; // Critical for getpar

  double *temp_par = (double *)mallocq((size_t)input.npar * sizeof(double));
  if (!temp_par)
  {
    puts("mallocq failed for temp_par");
    free(temp_idpar);
    fclose(lubak);
    return false;
  }
  double *temp_erp = (double *)mallocq((size_t)input.npar * sizeof(double));
  if (!temp_erp)
  {
    puts("mallocq failed for temp_erp");
    free(temp_idpar);
    free(temp_par);
    fclose(lubak);
    return false;
  }

  size_t parlbl_size = (LBLEN * (size_t)input.npar + 1);
  char *temp_parlbl = (char *)mallocq(parlbl_size);
  if (!temp_parlbl)
  {
    puts("mallocq failed for temp_parlbl");
    free(temp_idpar);
    free(temp_par);
    free(temp_erp);
    fclose(lubak);
    return false;
  }

  int original_npar_for_getpar = input.npar; // getpar can modify npar
  input.inpcor = getpar(lubak, lufit_for_getpar, &input.nfit, &original_npar_for_getpar,
                        temp_idpar, temp_par, temp_erp, temp_parlbl, LBLEN);
  input.npar = original_npar_for_getpar; // Update npar with value possibly changed by getpar

  // Copy data to CalFitInput's vectors
  // Adjust idpar_elem_count if npar changed in getpar, though it shouldn't affect the BCD definition part
  idpar_elem_count = (size_t)input.npar * input.ndbcd_from_options + input.ndbcd_from_options + 3;
  input.idpar_data.assign(temp_idpar, temp_idpar + idpar_elem_count);
  input.par_initial.assign(temp_par, temp_par + input.npar);
  input.erp_initial.assign(temp_erp, temp_erp + input.npar);
  parlbl_size = (LBLEN * (size_t)input.npar + 1); // Recalculate if npar changed
  input.parlbl_data_flat.assign(temp_parlbl, temp_parlbl + parlbl_size);

  free(temp_idpar);
  free(temp_par);
  free(temp_erp);
  free(temp_parlbl);

  // Handle getvar
  if (input.nfit > 0)
  {
    size_t nlsq_var_count = ((size_t)input.nfit * ((size_t)input.nfit + 1)) / 2;
    double *temp_var = (double *)mallocq(nlsq_var_count * sizeof(double));
    if (!temp_var)
    {
      puts("mallocq failed for temp_var");
      fclose(lubak);
      return false;
    }
    memset(temp_var, 0, nlsq_var_count * sizeof(double)); // Important for default case in getvar

    // idpar and erpar need to be passed as non-const pointers if getvar might modify them (unlikely for these)
    // However, .data() from const vector is const*. Use copies if modification is possible.
    // Assuming getvar doesn't modify idpar content pointed to, nor erpar content.
    // If it does, we need mutable copies.
    // For now, assume read-only access to idpar_data and erp_initial by getvar.
    input.inpcor = getvar(lubak, input.nfit, temp_var,
                          input.idpar_data.data(),  // Pass pointer to underlying data
                          input.erp_initial.data(), // Pass pointer to underlying data
                          input.inpcor);

    input.var_initial_from_getvar.assign(temp_var, temp_var + nlsq_var_count);
    free(temp_var);
  }
  else
  {
    input.var_initial_from_getvar.clear(); // No fit parameters, no variance matrix
  }

  // Read lines from the line file into input.lineData_raw
  FILE *lulin = fopen(linFile.c_str(), "r");
  if (!lulin)
  {
    printf("Error: Unable to open line file '%s'.\n", linFile.c_str());
    fclose(lubak);
    return false;
  }
  printf("Successfully opened line file '%s' for reading.\n", linFile.c_str());
  input.lineData_raw.clear();
  char lineBuffer[NDCARD]; // Assuming NDCARD is sufficient for line file lines too
  while (fgetstr(lineBuffer, NDCARD, lulin) > 0)
  {
    input.lineData_raw.push_back(std::string(lineBuffer));
  }
  printf("Read %zu lines from line file '%s' into memory.\n", input.lineData_raw.size(), linFile.c_str());

  fclose(lulin);
  fclose(lubak);

  return true;
}

// CalFitIO::writeOutput implementation will be complex and deferred after CalFit.run works
bool CalFitIO::writeOutput(const std::string &par_filepath_final,
                           const std::string &bak_filepath_original,
                           const std::string &var_filepath_final,
                           const CalFitOutput &output,
                           const CalFitInput &original_input)
{
  // This is a placeholder for Step 4 or later.
  // It will need to reconstruct the .par, .var files accurately.
  // For now, just create empty files or minimal output to show flow.
  FILE *lupar_final = fopen(par_filepath_final.c_str(), "w");
  if (lupar_final)
  {
    fprintf(lupar_final, "Title: %s", output.title_for_output.c_str());
    // ... more comprehensive writing ...
    fclose(lupar_final);
  }
  else
  {
    printf("Error: Unable to open final parameter file '%s' for writing.\n", par_filepath_final.c_str());
    return false;
  }

  FILE *luvar_final = fopen(var_filepath_final.c_str(), "w");
  if (luvar_final)
  {
    fprintf(luvar_final, "Title: %s", output.title_for_output.c_str());
    // ... more comprehensive writing, including putvar ...
    fclose(luvar_final);
  }
  else
  {
    printf("Error: Unable to open final variance file '%s' for writing.\n", var_filepath_final.c_str());
    return false;
  }

  // The bak_filepath_original is the original .par file, usually just kept as is or handled by filbak.
  // filbak in main.cpp handles the backup creation. This function might not need bak_filepath_original.

  printf("CalFitIO::writeOutput - Placeholder for writing to %s and %s\n",
         par_filepath_final.c_str(), var_filepath_final.c_str());

  return true;
}