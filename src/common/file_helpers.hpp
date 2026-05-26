/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

#ifndef FILE_HELPERS_HPP
#define FILE_HELPERS_HPP

#include <cstdio>

namespace file_helpers {

/* Open file for reading; search CWD then SPECNAME env path.
   Returns FILE* on success, nullptr if not found (no error thrown). */
FILE *open_input_optional(const char *fname);

/* Open file for reading; search CWD then SPECNAME env path.
   Throws IoError if not found. */
FILE *open_input(const char *fname);

/* Open file for writing (or other non-read mode).
   Throws IoError on failure. */
FILE *open_output(const char *fname, const char *mode = "w");

/* Parse file arguments from argc/argv (replaces filget).
   Same semantics: fills cfil[0..nfile-1] from argv extensions and cfil[nfile]
   with options string.  Returns nfile on success, 0 if no master found. */
int parse_file_args(int argc, char *argv[], int nfile,
                    char *cfil[], const char *cext[]);

/* Returns the base name (no extension) from a path, allocated with malloc().
   E.g. "mol.par" → "mol", "path/mol.toml" → "path/mol".
   Caller must free() the result. */
char *base_name(const char *path);

/* True if a file with the given path exists and is readable. */
bool file_exists(const char *path);

} // namespace file_helpers

#endif // FILE_HELPERS_HPP
