#ifndef OUTPUTSINK_HPP
#define OUTPUTSINK_HPP

#include <cstdarg>
#include <cstdio>
#include <string>
#include <vector>

/**
 * @brief Abstract output sink replacing raw FILE* in CalCat's public API.
 *
 * FileSink wraps a FILE* for CLI use. MemorySink buffers output in memory
 * for library/Python use, avoiding POSIX-only open_memstream/fmemopen.
 */
class OutputSink {
public:
    virtual ~OutputSink() = default;

    void printf(const char *fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        vprintf_impl(fmt, ap);
        va_end(ap);
    }

    // Write string without adding a newline (like fputs)
    virtual void puts(const char *s) = 0;

    // vprintf implementation — delegates from printf()
    virtual void vprintf_impl(const char *fmt, va_list ap) = 0;

    // Returns the underlying FILE* for legacy C interfaces (e.g. engine calls).
    // MemorySink returns a lazy tmpfile() bit-bucket; FileSink returns the real file.
    virtual FILE *file() const = 0;
};

// ---------------------------------------------------------------------------

class FileSink : public OutputSink {
    FILE *f_;
public:
    explicit FileSink(FILE *f) : f_(f) {}
    void puts(const char *s) override { fputs(s, f_); }
    void vprintf_impl(const char *fmt, va_list ap) override { vfprintf(f_, fmt, ap); }
    FILE *file() const override { return f_; }
};

// ---------------------------------------------------------------------------

class MemorySink : public OutputSink {
    std::string buf_;
    mutable FILE *fallback_; // lazy tmpfile for legacy FILE* callers

public:
    MemorySink() : fallback_(nullptr) {}
    ~MemorySink() override { if (fallback_) fclose(fallback_); }

    void puts(const char *s) override { buf_ += s; }

    void vprintf_impl(const char *fmt, va_list ap) override {
        va_list ap2;
        va_copy(ap2, ap);
        int n = vsnprintf(nullptr, 0, fmt, ap2);
        va_end(ap2);
        if (n > 0) {
            size_t pos = buf_.size();
            buf_.resize(pos + static_cast<size_t>(n) + 1);
            vsnprintf(&buf_[pos], static_cast<size_t>(n) + 1, fmt, ap);
            buf_.resize(pos + static_cast<size_t>(n));
        }
    }

    // Returns a tmpfile bit-bucket; engine writes here are not captured in buf_.
    FILE *file() const override {
        if (!fallback_) fallback_ = tmpfile();
        return fallback_;
    }

    // Splits buffered content on '\n', returns lines (without newline), clears buffer.
    std::vector<std::string> drain_lines() {
        std::vector<std::string> lines;
        size_t start = 0, end;
        while ((end = buf_.find('\n', start)) != std::string::npos) {
            lines.push_back(buf_.substr(start, end - start));
            start = end + 1;
        }
        if (start < buf_.size())
            lines.push_back(buf_.substr(start));
        buf_.clear();
        return lines;
    }

    const std::string &buffer() const { return buf_; }
};

#endif // OUTPUTSINK_HPP
