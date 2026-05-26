#ifndef TOML_IO_HPP
#define TOML_IO_HPP

/**
 * toml_io.hpp — TOML file format for SPFIT/SPCAT
 *
 * File roles:
 *   mol.toml       FitInput  — user-authored, read by spfit
 *   mol.var.toml   FitOutput — written by spfit, read by spcat
 *   mol.int.toml   CatInput extras (control + dipoles)
 *
 * These functions are linked into the CLI binaries (spfit/spcat) only.
 * They depend on toml++ (third_party/tomlplusplus/toml.hpp).
 */

#include "api/InputSchema.hpp"
#include "spfit/CalFit.hpp"
#include "spcat/CalCat.hpp"
#include <string>

/** Read mol.toml → FitInput.  Throws IoError / InputError on failure. */
FitInput  load_fit_input_toml(const std::string &path);

/** Read mol.var.toml + mol.int.toml → CatInput. */
CatInput  load_cat_input_toml(const std::string &var_path, const std::string &int_path);

/** Write CalFitOutput → mol.var.toml (merges engine_options + labels from src_fi). */
void save_fit_output_toml(const CalFitOutput &out, const FitInput &src_fi,
                          const std::string &path);

/** Write CalCatOutput → mol.cat.toml. */
void save_cat_output_toml(const CalCatOutput &out, const std::string &path);

#endif // TOML_IO_HPP
