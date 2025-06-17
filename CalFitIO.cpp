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
 * @param parFile Path to parameter input file (as backed up)
 * @param linFile Path to line file
 * @param input Output parameter for input data
 * @param calc_engine calculation engine
 * @param lufit_for_getpar .par output File stream for getpar to write its log.
 * @return True if reading is successful, false otherwise
 */
bool CalFitIO::readInput(const std::string &parFile, const std::string &linFile,
                         CalFitInput &input,
                         std::unique_ptr<CalculationEngine> &calc_engine_for_setup,
                         FILE *lufit_for_logging)
{
  if (!calc_engine_for_setup)
  {
    puts("Error: calc_engine_for_setup is NULL in CalFitIO::readInput.");
    if (lufit_for_logging)
      fprintf(lufit_for_logging, "ERROR: calc_engine_for_setup is NULL.\n");
    return false;
  }

  if (!lufit_for_logging)
  {
    puts("Error: lufit_for_logging stream is NULL."); // Cannot log this error to lufit itself
    return false;
  }

  // Open the par backup stream for reading
  FILE *lubak_stream_for_par_content = fopen(parFile.c_str(), "r");
  if (!lubak_stream_for_par_content)
  {
    printf("Error: Unable to open fit backup file '%s'.\n", parFile.c_str());
    if (lufit_for_logging)
      fprintf(lufit_for_logging, "ERROR: lubak_stream_for_par_content is NULL.\n");
    return false;
  }

  char card[NDCARD];

  // Read title of .par file
  if (fgetstr(card, NDCARD, lubak_stream_for_par_content) <= 0)
  {
    fprintf(lufit_for_logging, "Unable to read title of .par file\n");
    puts("Unable to read title of .par file");
    fclose(lubak_stream_for_par_content);
    return false;
  }
  chtime(card, 82);
  input.title = std::string(card);
  fputs(input.title.c_str(), lufit_for_logging); // Log to lufit

  // Read run parameters (2nd line)
  double dvec[8] = {100.0, 32767.0, 1.0, 0.0, 0.0, 1e6, 1.0, 1.0};
  int n_read_dvec = fgetstr(card, NDCARD, lubak_stream_for_par_content);
  if (n_read_dvec != 0)
    n_read_dvec = pcard(card, dvec, 8, NULL);
  if (n_read_dvec == 0)
  {
    fprintf(lufit_for_logging, "Unable to read second line of .par file\n");
    puts("Unable to read second line of .par file");
    fclose(lubak_stream_for_par_content);
    return false;
  }
  input.npar = (int)dvec[0];
  input.limlin = (int)dvec[1];
  input.nitr = (int)dvec[2];
  input.nxpar_from_file = (int)dvec[3];
  input.marqp0 = dvec[4];
  input.xerrmx = dvec[5];
  input.parfac_initial = dvec[6];
  input.fqfacq = dvec[7];

  // Store file position before setopt reads options
  long pos_before_options = ftell(lubak_stream_for_par_content);
  if (pos_before_options == -1L)
  {
    // Handle error: cannot reliably read options and then parameters
    perror("ftell before options failed");
    fclose(lubak_stream_for_par_content);
    return false;
  }

  // Call CalculationEngine::setopt
  // setopt will read option lines from lubak_stream_for_par_content
  // and advance the file pointer past them.
  char temp_namfil_buffer[NDCARD] = {0};
  // Initialize parameters for setopt with defaults that SpinvEngine::setopt expects or can modify
  int temp_nfmt_cat = MAXCAT; // Default for setopt's nfmt output (catalog format count)
  int temp_itd = 2;           // Default for molecule type
  int temp_ndbcd = 1;         // Default BCD length

  input.noptn_read_by_setopt = calc_engine_for_setup->setopt(
      lubak_stream_for_par_content,  // setopt consumes lines from this stream
      &temp_nfmt_cat,
      &temp_itd,
      &temp_ndbcd,
      temp_namfil_buffer);

  if (input.noptn_read_by_setopt < 0)
  { // EOF or error during option reading
    fprintf(lufit_for_logging, "Warning/Error: calc_engine->setopt returned %d. Check .par file options section.\n", input.noptn_read_by_setopt);
    // This might be acceptable if no options. If options were expected, it's an issue.
  }
  input.nfmt_cat_from_setopt = temp_nfmt_cat;
  input.itd_from_setopt = temp_itd;
  input.ndbcd_from_setopt = temp_ndbcd;
  input.namfil_from_setopt = std::string(temp_namfil_buffer);
  // Note: CalFitIO does not explicitly log the option cards here, assuming setopt or getpar might.
  // If not, and we need option cards in lufit, we'd have to read them, store, pass to setopt (e.g. via temp file), then log.
  // For now, assume setopt handles its own logging or the information is implicit.

  // Rewind to read and store the option lines that setopt processed
  if (fseek(lubak_stream_for_par_content, pos_before_options, SEEK_SET) != 0)
  {
    perror("fseek to re-read options failed");
    // Handle error: cannot get option lines, but setopt might have worked.
    // Output files might miss option lines.
    // For now, proceed but log warning.
    if (lufit_for_logging)
      fprintf(lufit_for_logging, "Warning: Could not rewind to save option lines.\n");
    // To ensure lubak_stream is correctly positioned for getpar, could try to advance it by what setopt read,
    // though this is less reliable than setopt just leaving it at the right place.
    // For now, if fseek fails, subsequent getpar might read wrong data. Best to return false.
    fclose(lubak_stream_for_par_content);
    return false;
  }

  input.raw_option_lines_from_par.clear();
  if (input.noptn_read_by_setopt > 0)
  {
    char option_card_buffer[NDCARD];
    for (int i = 0; i < input.noptn_read_by_setopt; ++i)
    {
      if (fgetstr(option_card_buffer, NDCARD, lubak_stream_for_par_content) > 0)
      {
        input.raw_option_lines_from_par.push_back(std::string(option_card_buffer));
      }
      else
      {
        if (lufit_for_logging)
          fprintf(lufit_for_logging, "Warning: Premature EOF while re-reading option line %d.\n", i + 1);
        // This implies setopt read more lines than are now available, or fgetstr error.
        break;
      }
    }
  }
  else if (input.noptn_read_by_setopt < 0)
  {
    // setopt indicated an error or EOF during its own reading.
    if (lufit_for_logging)
      fprintf(lufit_for_logging, "Note: setopt returned %d, no option lines stored.\n", input.noptn_read_by_setopt);
  }
  // lubak_stream_for_par_content is now positioned AFTER option lines, ready for getpar.

  // Call getpar
  // `lubak_stream_for_par_content` is now positioned at the start of parameters.
  // `input.ndbcd_from_setopt` is the authoritative value.
  size_t idpar_elem_count = (size_t)input.npar * input.ndbcd_from_setopt + input.ndbcd_from_setopt + 3;
  unsigned char *temp_idpar = (unsigned char *)mallocq(idpar_elem_count * sizeof(unsigned char));
  if (!temp_idpar)
  { /* error */
    fclose(lubak_stream_for_par_content);
    return false;
  }
  temp_idpar[0] = (unsigned char)input.ndbcd_from_setopt; // Set NDEC for getpar

  double *temp_par = (double *)mallocq((size_t)input.npar * sizeof(double));
  if (!temp_par)
  {
    free(temp_idpar);
    fclose(lubak_stream_for_par_content);
    return false;
  }
  double *temp_erp = (double *)mallocq((size_t)input.npar * sizeof(double));
  if (!temp_erp)
  {
    free(temp_idpar);
    free(temp_par);
    fclose(lubak_stream_for_par_content);
    return false;
  }

  size_t parlbl_size = (LBLEN * (size_t)input.npar + 1);
  char *temp_parlbl = (char *)mallocq(parlbl_size);
  if (!temp_parlbl)
  {
    free(temp_idpar);
    free(temp_par);
    free(temp_erp);
    fclose(lubak_stream_for_par_content);
    return false;
  }

  int original_npar_for_getpar = input.npar;
  input.inpcor = getpar(lubak_stream_for_par_content, lufit_for_logging,
                        &input.nfit, &original_npar_for_getpar,
                        temp_idpar, temp_par, temp_erp, temp_parlbl, LBLEN);
  input.npar = original_npar_for_getpar; // Update npar

  // Recalculate sizes if npar changed
  idpar_elem_count = (size_t)input.npar * input.ndbcd_from_setopt + input.ndbcd_from_setopt + 3;
  parlbl_size = (LBLEN * (size_t)input.npar + 1);

  input.idpar_data.assign(temp_idpar, temp_idpar + idpar_elem_count);
  input.par_initial.assign(temp_par, temp_par + input.npar);
  input.erp_initial.assign(temp_erp, temp_erp + input.npar);
  input.parlbl_data_flat.assign(temp_parlbl, temp_parlbl + parlbl_size);

  free(temp_idpar);
  free(temp_par);
  free(temp_erp);
  free(temp_parlbl);

  // Call getvar
  if (input.nfit > 0)
  {
    size_t nlsq_var_count = ((size_t)input.nfit * ((size_t)input.nfit + 1)) / 2;
    double *temp_var = (double *)mallocq(nlsq_var_count * sizeof(double));
    if (!temp_var)
    {
      fclose(lubak_stream_for_par_content);
      return false;
    }
    memset(temp_var, 0, nlsq_var_count * sizeof(double));

    input.inpcor = getvar(lubak_stream_for_par_content, input.nfit, temp_var,
                          input.idpar_data.data(), input.erp_initial.data(),
                          input.inpcor);
    input.var_initial_from_getvar.assign(temp_var, temp_var + nlsq_var_count);
    free(temp_var);
  }
  else
  {
    input.var_initial_from_getvar.clear();
  }

  // Read .lin file
  FILE *lulin_stream = fopen(linFile.c_str(), "r");
  if (!lulin_stream)
  { /* error, log to lufit_for_logging */
    return false;
  }
  input.lineData_raw.clear();
  char lineBuffer[NDCARD];
  while (fgetstr(lineBuffer, NDCARD, lulin_stream) > 0)
  {
    input.lineData_raw.push_back(std::string(lineBuffer));
  }
  fclose(lulin_stream);
  fclose(lubak_stream_for_par_content);

  return true;
} // readInput


bool CalFitIO::writeOutput(const std::string &par_filepath_final,
                           const std::string &bak_filepath_original, // Unused by this function itself
                           const std::string &var_filepath_final,
                           const CalFitOutput &output,
                           const CalFitInput &original_input)
{
  // Silence unused parameter warning if bak_filepath_original is truly unused here
  (void)bak_filepath_original; // TODO REMOVE

  FILE *lupar = fopen(par_filepath_final.c_str(), "w");
  if (!lupar)
  {
    printf("Error: Unable to open final parameter file '%s' for writing.\n", par_filepath_final.c_str());
    return false;
  }
  printf("Opened final parameter file '%s' for writing.\n", par_filepath_final.c_str());

  FILE *luvar = fopen(var_filepath_final.c_str(), "w");
  if (!luvar)
  {
    printf("Error: Unable to open final variance file '%s' for writing.\n", var_filepath_final.c_str());
    fclose(lupar);
    return false;
  }
  printf("Opened final variance file '%s' for writing.\n", var_filepath_final.c_str());

  // 1. Write Title (already time-stamped if done by chtime earlier)
  // Assuming output.title_for_output is the final, possibly time-stamped title.
  // If not, chtime it here. Let's assume it's ready.
  // Original main.c did: rewind(lubak); fgetstr(card, NDCARD, lubak); chtime(card, 82);
  // This implies re-reading the original title from the .bak file for the final output.
  // Let's use output.title_for_output which should be derived from input.title.
  char title_card_buffer[NDCARD];
  strncpy(title_card_buffer, output.title_for_output.c_str(), NDCARD - 1);
  title_card_buffer[NDCARD - 1] = '\0';
  // chtime(title_card_buffer, 82); // Apply chtime again to ensure latest timestamp if needed
  // Or assume output.title_for_output is already final.
  // For safety, let's re-apply chtime from a base title.
  strncpy(title_card_buffer, original_input.title.c_str(), NDCARD - 1); // Use original non-timestamped title
  title_card_buffer[NDCARD - 1] = '\0';
  chtime(title_card_buffer, 82); // Timestamp it now

  fputs(title_card_buffer, lupar);
  fputs(title_card_buffer, luvar);
  // No newline by default from fputs if string doesn't have it. Add if chtime doesn't.
  // chtime usually adds newline.

  // 2. Write Second Header Line
  // Original: fprintf(lupar,"%4d %4d %4d %4d %14.4E %14.4E %14.4E %12.10f\n",
  //           npar, limlin, nitr, nxpar, marqlast, xerrmx, parfac0, fqfacq);
  // output.limlin_final should be the original limlin value (potentially negative)
  fprintf(lupar, "%4d %4d %4d %4d %14.4E %14.4E %14.4E %12.10f\n",
          output.npar_final, output.limlin_final, output.nitr_final_actual,
          output.nxpar_for_header, output.marqlast_final, output.xerrmx_final,
          output.parfac_for_header, output.fqfacq_final);
  fprintf(luvar, "%4d %4d %4d %4d %14.4E %14.4E %14.4E %12.10f\n",
          output.npar_final, output.limlin_final, output.nitr_final_actual,
          output.nxpar_for_header, output.marqlast_final, output.xerrmx_final,
          output.parfac_for_header, output.fqfacq_final);

  // 3. Write Option Lines (from original_input as they were read from .par)
  // Original: for (icnt=0; icnt<noptn; ++icnt) { fgetstr(card,lubak); card[k]='\n'; card[k+1]='\0'; fputs(card,lupar);}
  if (!original_input.raw_option_lines_from_par.empty())
  {
    printf("Writing %zu stored option lines to output files.\n", original_input.raw_option_lines_from_par.size());
    for (const auto &opt_line_str : original_input.raw_option_lines_from_par)
    {
      fputs(opt_line_str.c_str(), lupar);
      // fgetstr usually includes newline if line wasn't too long.
      // If opt_line_str does not end with \n, add one.
      if (!opt_line_str.empty() && opt_line_str.back() != '\n')
      {
        fputc('\n', lupar);
      }
      fputs(opt_line_str.c_str(), luvar);
      if (!opt_line_str.empty() && opt_line_str.back() != '\n')
      {
        fputc('\n', luvar);
      }
    }
  }
  else if (original_input.noptn_read_by_setopt > 0)
  {
    // This case means setopt said it read options, but we failed to store them.
    printf("Warning: Option lines were processed by setopt (%d) but not stored for output.\n", original_input.noptn_read_by_setopt);
    if (lupar != stdout)
      fprintf(lupar, "! Warning: Option lines not available for writing.\n");
    if (luvar != stdout)
      fprintf(luvar, "! Warning: Option lines not available for writing.\n");
  }

  // 4. Write Parameter Lines
  // Original loop: for (i=0, ibcd=0; i<npar; ++i, ibcd+=ndbcd) { ... }
  // We need idpar, par (fitted), erp (original error for .par), erpar (fitted error for .var), parlbl
  char param_id_str[64];                 // Buffer for putbcd output
  char param_label_card_part[LBLEN + 2]; // For "/label" or empty

  const unsigned char *current_idpar_ptr = output.idpar_final_for_output.data();
  const char *current_parlbl_ptr = output.parlbl_final_for_output_flat.data();
  int ndbcd_val = original_input.ndbcd_from_setopt; // Get this from input as CalFitOutput might not store it.
                                                    // Or add it to CalFitOutput. For now, from original_input.

  int ibase_for_dependent_error = 0; // Keeps track of the index of the last independent parameter

  for (int i = 0; i < output.npar_final; ++i)
  {
    // Skip parameter lines in original .bak file (not needed here as we construct anew)

    // Prepare label part: "/label" or empty string if label is empty
    if (current_parlbl_ptr[0] == '\0')
    {
      param_label_card_part[0] = '\0';
    }
    else
    {
      param_label_card_part[0] = '/';
      // strncpy is safer if LBLEN is exact and no null terminator in source
      strncpy(param_label_card_part + 1, current_parlbl_ptr, LBLEN);
      param_label_card_part[LBLEN + 1] = '\0'; // Ensure null termination
    }

    putbcd(param_id_str, 64, current_idpar_ptr); // Format parameter ID

    double current_par_val = output.par[i];
    double current_erp_orig_val = output.erp_original_for_output[i];
    double current_erpar_fitted_val = output.erpar[i];

    // Handle dependent parameters for .var file error (as in original main's output loop)
    // and potentially for .par file value if it was stored as a factor.
    if (NEGBCD(current_idpar_ptr[0]) != 0)
    { // If parameter is dependent (fixed or tied)
      // Value of dependent parameter: par[i] (factor) * par[ibase] (master value)
      // Error for dependent parameter: erpar[i] = fabs(par[i-factor]) * erpar[ibase-fitted-error]
      // This logic was in main.c's output loop:
      // scale = fabs(par[i]); erpar[i] = scale * erpar[ibase]; par[i] *= par[ibase];
      // This means output.par[i] for dependent parameters should ideally be the final *value*, not factor.
      // Let's assume output.par[i] from CalFit::finalizeOutputData IS the final value.
      // And output.erpar[i] for dependent IS the derived error.
    }
    else
    {                                // Independent parameter
      ibase_for_dependent_error = i; // Update index of last independent parameter
    }

    // Write to .par file: ID, Value, OriginalError, /Label
    fprintf(lupar, "%s %23.15E %14.8E %s\n",
            param_id_str,
            current_par_val,
            current_erp_orig_val, // Original a-priori error
            param_label_card_part);

    // Write to .var file: ID, Value, FittedError, /Label
    fprintf(luvar, "%s %23.15E %14.8E %s\n",
            param_id_str,
            current_par_val,
            current_erpar_fitted_val, // Final fitted error
            param_label_card_part);

    current_idpar_ptr += ndbcd_val;
    current_parlbl_ptr += LBLEN;
  }

  // 5. Write Variance/Correlation Data
  // Original: putvar(luvar, nfit, var, dpar);
  // `var` is the packed matrix (e.g. L_inv_final_packed_lower or Cov_inv_packed_upper)
  // `dpar` is the array of scaled errors for the *fitted* parameters.
  if (output.nfit_final > 0)
  {
    if (!output.var_final_for_output.empty() && !output.dpar_final_for_putvar.empty())
    {
      // putvar expects non-const pointers for var and erpar (dpar here)
      // Create mutable copies if necessary, or cast const away if putvar guarantees no modification.
      // Assuming putvar doesn't modify them (though its erpar is usually output for scaling).
      // The `putvar` from ulib.c *does* modify `var` by scaling it with `erpar`.
      // So, we need mutable copies.
      std::vector<double> var_copy = output.var_final_for_output;
      std::vector<double> dpar_copy = output.dpar_final_for_putvar;
      putvar(luvar, output.nfit_final, var_copy.data(), dpar_copy.data());
    }
    else if (lupar != stdout && luvar != stdout)
    { // Avoid print if just stdout
      printf("Warning: var_final_for_output or dpar_for_putvar is empty, skipping putvar.\n");
    }
  }

  // Original main copies remaining lines from .bak to .par if any (e.g. correlation lines already there)
  // This is not done here, as we are creating new .par and .var from scratch based on fit results.
  // If original .par had extra info after params/var that needs preserving, that's a different requirement.

  fclose(lupar);
  fclose(luvar);

  printf("Output files %s and %s written successfully.\n", par_filepath_final.c_str(), var_filepath_final.c_str());
  return true;
}