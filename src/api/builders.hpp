#ifndef BUILDERS_HPP
#define BUILDERS_HPP

#include "api/InputSchema.hpp"
#include "spfit/CalFit.hpp"
#include "spcat/CalCat.hpp"
#include "engine/CalculationEngine.hpp"
#include "common/Logger.hpp"

/**
 * Build a CalFitInput from a typed FitInput struct.
 *
 * Calls eng.apply_options() to configure the engine, BCD-packs parameter IDs,
 * converts absolute parameter values to ratios for dependent parameters, and
 * formats line records into the legacy string format that CalFit::linein expects.
 *
 * Throws ValidationError / InputError on schema violations.
 */
CalFitInput build_fit_input(const FitInput &fi, CalculationEngine &eng,
                            Logger &logger = Logger::defaultLogger());

/**
 * Build a CalCatInput from a typed CatInput struct.
 *
 * Calls eng.apply_options(), BCD-packs dipole and parameter IDs, and fills
 * the packed upper-triangular variance matrix.
 *
 * Throws ValidationError / InputError on schema violations.
 */
CalCatInput build_cat_input(const CatInput &ci, CalculationEngine &eng,
                            Logger &logger = Logger::defaultLogger());

#endif // BUILDERS_HPP
