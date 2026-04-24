#include "file_helpers.hpp"
#include "CalError.hpp"
#include <cstdlib>
#include <cstring>
#include <string>

namespace file_helpers {

FILE *open_input_optional(const char *fname)
{
    FILE *h = fopen(fname, "r");
    if (h == nullptr) {
        const char *path = getenv("SPECNAME");
        if (path != nullptr) {
            std::string full = std::string(path) + fname;
            h = fopen(full.c_str(), "r");
        } else {
            puts("path variable SPECNAME can be set to path for names");
        }
    }
    return h;
}

FILE *open_output(const char *fname, const char *mode)
{
    FILE *h = fopen(fname, mode);
    if (h == nullptr)
        throw IoError(std::string("Trouble opening ") + fname);
    return h;
}

int parse_file_args(int argc, char *argv[], int nfile,
                    char *cfil[], const char *cext[])
{
#define NSTR 82
    char str[NSTR];
    char *master = nullptr, *opt = nullptr, *parg, *ext, *pstr;
    int iext, k, len = 0, iarg, knt, query;

    for (iext = 0; iext <= nfile; ++iext) cfil[iext] = nullptr;
    knt = argc - 1; query = 0;
    if (knt <= 0) { query = 1; knt = nfile + 1; }

    for (iarg = 1; iarg <= knt; ++iarg) {
        if (query) {
            parg = str;
            puts(" Enter file name ");
            /* use fgetstr semantics inline */
            if (fgets(parg, NSTR, stdin) == nullptr) break;
            int n = 0, ki = 0;
            for (char *p = parg; *p && *p != '\n'; ++p, ++ki)
                if (*p != ' ' && (unsigned char)*p >= 0x20) n = ki + 1;
            parg[n] = '\0';
            if (n <= 0) continue;
        } else {
            parg = argv[iarg];
            if (parg == nullptr) continue;
        }
        ext = strrchr(parg, '.');
        if (ext != nullptr) {
            ++ext;
            for (iext = 0; iext < nfile; ++iext) {
                if (strcmp(cext[iext], ext) == 0) {
                    k = (int)strlen(parg) + 1;
                    pstr = (char *)malloc((size_t)k);
                    if (pstr) memcpy(pstr, parg, (size_t)k);
                    cfil[iext] = pstr;
                    break;
                }
            }
        }
        if (parg[0] == '-') {
            if (opt) continue;
            k = (int)strlen(parg);
            opt = (char *)malloc((size_t)k);
            if (!opt) break;
            --k;
            memcpy(opt, parg + 1, (size_t)k);
            opt[k] = '\0'; cfil[nfile] = opt;
            continue;
        }
        if (!master) {
            k = (int)strlen(parg) + 1;
            if (ext) k -= (int)strlen(ext) + 1;
            if (k <= 0) continue;
            master = (char *)malloc((size_t)k);
            if (!master) break;
            len = k; --k;
            memcpy(master, parg, (size_t)k);
            master[k] = '.';
        }
    }
    if (!master) return 0;
    for (iext = 0; iext < nfile; ++iext) {
        pstr = cfil[iext];
        if (!pstr) {
            k = (int)strlen(cext[iext]) + 1;
            pstr = (char *)malloc((size_t)(len + k));
            if (pstr) {
                memcpy(pstr, master, (size_t)len);
                memcpy(pstr + len, cext[iext], (size_t)k);
                cfil[iext] = pstr;
            }
        }
        if (pstr) puts(pstr);
    }
    free(master);
    return nfile;
#undef NSTR
}

} // namespace file_helpers
