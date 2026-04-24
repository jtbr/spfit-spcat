/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

#ifndef CALERROR_HPP
#define CALERROR_HPP

#include <stdexcept>
#include <string>

enum class CalErrorCode {
    Unknown = 0,
    // I/O
    FileOpenFailed,
    FileReadFailed,
    FileWriteFailed,
    // Malformed input data
    MalformedInput,
    // API-boundary structural violations
    NullInput,
    SizeMismatch,
    InvalidParameter,
    // Numerical failures
    DiagonalizationFailed,
    ScratchFileIo,
    WorkingVectorTooShort,
    AllocationFailed,
    // Physics/symmetry violations caught at runtime
    SpinDimensioning,
    LDoublets,
};

class CalError : public std::runtime_error {
public:
    CalError(CalErrorCode code, const std::string &msg)
        : std::runtime_error(msg), m_code(code) {}
    CalErrorCode code() const { return m_code; }
private:
    CalErrorCode m_code;
};

/* File-open / file-read / file-write failures */
class IoError : public CalError {
public:
    explicit IoError(const std::string &msg,
                     CalErrorCode code = CalErrorCode::FileOpenFailed)
        : CalError(code, msg) {}
};

/* Malformed .par / .lin / .int / .var content */
class InputError : public CalError {
public:
    explicit InputError(const std::string &msg,
                        CalErrorCode code = CalErrorCode::MalformedInput)
        : CalError(code, msg) {}
};

/* API-boundary structural violations (wrong sizes, null pointers, etc.) */
class ValidationError : public CalError {
public:
    explicit ValidationError(const std::string &msg,
                             CalErrorCode code = CalErrorCode::NullInput)
        : CalError(code, msg) {}
};

/* Numerical failures: diagonalization, scratch I/O, work-vector bounds */
class NumericError : public CalError {
public:
    explicit NumericError(const std::string &msg,
                          CalErrorCode code = CalErrorCode::DiagonalizationFailed)
        : CalError(code, msg) {}
};

/* Allocate n bytes; throw NumericError on failure (C++ replacement for mallocq). */
#include <cstdlib>
inline void *calalloc(std::size_t n) {
    void *p = std::malloc(n);
    if (p == nullptr && n > 0)
        throw NumericError("memory allocation error", CalErrorCode::AllocationFailed);
    return p;
}

#endif // CALERROR_HPP
