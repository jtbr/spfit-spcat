/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <cstdio>
#include <string>

enum class LogLevel { debug, info, warn, error };

/* Abstract logger interface. Callers (and Python bindings) can install a custom
   sink; the default writes info/warn/error to stdout/stderr and suppresses debug. */
class Logger {
public:
    virtual ~Logger() = default;
    virtual void log(LogLevel level, const std::string &msg) = 0;

    /* Plain-string overloads */
    void debug(const std::string &msg) { log(LogLevel::debug, msg); }
    void info (const std::string &msg) { log(LogLevel::info,  msg); }
    void warn (const std::string &msg) { log(LogLevel::warn,  msg); }
    void error(const std::string &msg) { log(LogLevel::error, msg); }

    /* Format overloads — require at least one format argument so plain-string
       calls (no args) always go to the std::string overload above. */
    template<typename A, typename... R>
    void debug(const char *fmt, A a, R... r) { log(LogLevel::debug, format(fmt, a, r...)); }
    template<typename A, typename... R>
    void info (const char *fmt, A a, R... r) { log(LogLevel::info,  format(fmt, a, r...)); }
    template<typename A, typename... R>
    void warn (const char *fmt, A a, R... r) { log(LogLevel::warn,  format(fmt, a, r...)); }
    template<typename A, typename... R>
    void error(const char *fmt, A a, R... r) { log(LogLevel::error, format(fmt, a, r...)); }

    /* Returns the process-lifetime default logger (writes to stdout/stderr). */
    static Logger &defaultLogger();

private:
    template<typename A, typename... R>
    static std::string format(const char *fmt, A a, R... r) {
        int n = std::snprintf(nullptr, 0, fmt, a, r...);
        std::string s(static_cast<std::size_t>(n), '\0');
        std::snprintf(s.data(), static_cast<std::size_t>(n) + 1, fmt, a, r...);
        return s;
    }
};

/* Default implementation: info→stdout, warn/error→stderr, debug suppressed. */
class StdioLogger : public Logger {
public:
    void log(LogLevel level, const std::string &msg) override;
};

#endif // LOGGER_HPP
