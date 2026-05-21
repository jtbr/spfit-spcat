#ifndef LEGACY_PARSER_HPP
#define LEGACY_PARSER_HPP

#include "api/InputSchema.hpp"
#include <string>

/**
 * Parse a SPFIT .par + .lin file pair into a FitInput struct.
 *
 * Option cards are decoded into EngineOptions.spinv (or .dpi if kind == Dpi).
 * Parameter lines are decoded into FitInput.parameters (absolute values).
 * .lin lines are stored verbatim in FitInput.raw_lines for build_fit_input.
 * Variance is left empty (build_fit_input will build the default diagonal).
 *
 * Does NOT configure the engine — call build_fit_input(fi, engine) afterwards,
 * which calls apply_options and encodes parameters into CalFitInput.
 *
 * Throws InputError / IoError on file or parse failures.
 */
FitInput parse_fit_files(const std::string &parFile,
                         const std::string &linFile,
                         EngineOptions::Kind kind = EngineOptions::Kind::Spinv);

/**
 * Parse a SPCAT .var + .int file pair into a CatInput struct.
 *
 * Same conventions as parse_fit_files; variance is read from the .var file.
 *
 * Throws InputError / IoError on file or parse failures.
 */
CatInput parse_cat_files(const std::string &varFile,
                         const std::string &intFile,
                         EngineOptions::Kind kind = EngineOptions::Kind::Spinv);

#endif // LEGACY_PARSER_HPP
