/**
 * cat_example.cpp — Programmatic CalCat usage with in-memory output.
 *
 * Reads a SPCAT-format molecule, generates the spectroscopic catalog using
 * MemorySinks (no output files written), sorts the captured catalog lines by
 * frequency, and prints the first N lines plus the 300 K partition function.
 *
 * The .int file comes from <int_base>.int and the .var file from <var_base>.var.
 * The .var file is produced by spfit; point var_base at a pre-computed output
 * directory if spfit has not been run yet.
 *
 * Build:  cmake -DBUILD_EXAMPLES=ON .. && make cat_example
 * Usage:  cat_example <int_base> <var_base> [max_lines]
 * Example (from project root, using the v2008 pre-computed .var):
 *   cat_example  spfit_spcat_test_suite/diatomic_molecules/co_4/co_4  \
 *                spfit_spcat_test_suite/diatomic_molecules/co_4/v2008_results/co_4
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

#include "engine/SpinvEngine.hpp"
#include "spcat/CalCat.hpp"
#include "spcat/CalCatIO.hpp"
#include "spcat/OutputSink.hpp"
#include "common/CalError.hpp"


int main(int argc, char *argv[])
{
    if (argc < 3) {
        fprintf(stderr,
                "Usage: %s <int_base> <var_base> [max_lines]\n"
                "  int_base  — base path supplying <int_base>.int\n"
                "  var_base  — base path supplying <var_base>.var (spfit output)\n"
                "  max_lines — catalog lines to print (default 10)\n\n"
                "Example:\n"
                "  %s spfit_spcat_test_suite/diatomic_molecules/co_4/co_4 \\\n"
                "     spfit_spcat_test_suite/diatomic_molecules/co_4/v2008_results/co_4\n",
                argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    const std::string int_file = std::string(argv[1]) + ".int";
    const std::string var_file = std::string(argv[2]) + ".var";
    const int max_lines = (argc >= 4) ? std::atoi(argv[3]) : 10;

    // Diagnostic output (equivalent to .out file) goes to stdout.
    FileSink out_sink(stdout);

    std::unique_ptr<CalculationEngine> engine = std::make_unique<SpinvEngine>();

    CalCatInput input;
    if (!CalCatIO::readInput(int_file, var_file, input, engine, &out_sink)) {
        fprintf(stderr, "Error: could not read input files.\n");
        return EXIT_FAILURE;
    }

    // Capture .cat, .egy, .str output in memory instead of writing files.
    MemorySink cat_sink, egy_sink, str_sink;

    CalCat calCat(engine, &out_sink, &cat_sink, &egy_sink, &str_sink);
    CalCatOutput output;
    try {
        calCat.run(input, output);
    } catch (const CalError &e) {
        fprintf(stderr, "Catalog generation failed: %s\n", e.what());
        return EXIT_FAILURE;
    }

    // Drain in-memory sinks into the output struct.
    output.cat_lines = cat_sink.drain_lines();
    output.egy_lines = egy_sink.drain_lines();
    output.str_lines = str_sink.drain_lines();

    // Sort catalog lines by frequency (mimics the file-based sortn step).
    output.sort_cat_lines();

    // --- Programmatic access to results ---
    printf("\n=== CalCat results ===\n");
    printf("  Catalog lines generated : %ld\n", output.nline);

    int show = std::min(max_lines, (int)output.cat_lines.size());
    printf("  First %d catalog line(s):\n", show);
    for (int i = 0; i < show; ++i)
        printf("    %s\n", output.cat_lines[i].c_str());

    // Partition function at 300 K (closest temperature in the table).
    double q300 = 0.0, best_diff = 1e30;
    for (int i = 0; i < output.ntemp; ++i) {
        double diff = std::fabs(output.temp[i] - 300.0);
        if (diff < best_diff) { best_diff = diff; q300 = output.qsum[i]; }
    }
    printf("  Q(300 K)                : %.4f  (log10 = %.4f)\n",
           q300, (q300 > 0) ? std::log10(q300) : -100.0);

    return EXIT_SUCCESS;
}
