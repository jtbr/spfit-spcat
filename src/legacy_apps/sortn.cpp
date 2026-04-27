#include <stdio.h>
#include <stdlib.h>
#include "common/CalError.hpp"
#include "common/file_helpers.hpp"
#include "common/SigintFlag.hpp"
extern "C" int sortn(char *impname, char *outname, int dokey);

int main(int argc, char *argv[])
{
  try {
    FILE *finp;
    char *outnam;
    if (argc < 2) {
      puts("sortn infile [outfile]");
      exit(1);
    }
    finp = fopen(argv[1], "r");
    if (finp == NULL) {
      puts("no input file");
      exit(1);
    }
    if (argc >= 3) {
      outnam = argv[2];
    } else {
      outnam = argv[1];
    }
    return sortn(argv[1], outnam, 1);
  } catch (const CalError &e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return EXIT_FAILURE;
  } catch (const std::exception &e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return EXIT_FAILURE;
  }
}
