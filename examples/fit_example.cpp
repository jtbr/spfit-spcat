/**
 * fit_example.cpp — Programmatic CalFit usage demonstration.
 *
 * Reads a SPFIT-format molecule (base.par + base.lin), runs the parameter
 * fitting in-memory, and prints the fitted parameters and RMS error.
 * No output files are written — this exercises the pure library path.
 *
 * Build:  cmake -DBUILD_EXAMPLES=ON .. && make fit_example
 * Usage:  fit_example <base_path>
 * Example (from project root):
 *   fit_example spfit_spcat_test_suite/diatomic_molecules/co_4/co_4
 */

#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>

#include "engine/SpinvEngine.hpp"
#include "spfit/CalFit.hpp"
#include "spfit/CalFitIO.hpp"
#include "common/CalError.hpp"
#include "common/Logger.hpp"

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <base_path>\n", argv[0]);
        fprintf(stderr, "  e.g. %s spfit_spcat_test_suite/diatomic_molecules/co_4/co_4\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    const std::string base(argv[1]);
    const std::string par_file = base + ".par";
    const std::string lin_file = base + ".lin";

    // CalFit writes a fit log to a FILE*. We use stdout here so the normal
    // SPFIT output is visible; in a non-interactive context pass a tmpfile().
    FILE *fit_log = stdout;

    std::unique_ptr<CalculationEngine> engine = std::make_unique<SpinvEngine>();

    CalFitInput input;
    if (!CalFitIO::readInput(par_file, lin_file, input, engine, fit_log)) {
        fprintf(stderr, "Error: could not read input files '%s'.\n", base.c_str());
        return EXIT_FAILURE;
    }

    CalFit calFit(engine, fit_log); // engine ownership transferred here
    CalFitOutput output;
    try {
        calFit.run(input, output);
    } catch (const CalError &e) {
        fprintf(stderr, "Fit failed: %s\n", e.what());
        return EXIT_FAILURE;
    }

    // --- Programmatic access to results ---
    printf("\n=== CalFit results ===\n");
    printf("  RMS error  : %.6f\n", output.xsqbest);
    printf("  Iterations : %d\n",   output.itr);
    printf("  Parameters :\n");
    for (size_t i = 0; i < output.par.size(); ++i) {
        printf("    [%2zu]  %21.13E  +/-  %12.5E\n",
               i + 1, output.par[i], output.erpar[i]);
    }

    return EXIT_SUCCESS;
}
