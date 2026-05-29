#pragma once
// MSVC does not provide fmemopen (POSIX). This shim opens a temporary file,
// writes the buffer, and rewinds — identical semantics for our read-only use.
#ifdef _MSC_VER
#include <cstdio>
static inline FILE* fmemopen(const void* buf, size_t size, const char* /*mode*/) {
    FILE* f = tmpfile();
    if (!f) return nullptr;
    if (size > 0) {
        fwrite(buf, 1, size, f);
        rewind(f);
    }
    return f;
}
#endif
