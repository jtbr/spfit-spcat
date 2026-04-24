/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

enum class LogLevel { debug, info, warn, error };

/* Abstract logger interface. Callers (and Python bindings) can install a custom
   sink; the default writes info/warn/error to stdout/stderr and suppresses debug. */
class Logger {
public:
    virtual ~Logger() = default;
    virtual void log(LogLevel level, const std::string &msg) = 0;

    void debug(const std::string &msg) { log(LogLevel::debug, msg); }
    void info (const std::string &msg) { log(LogLevel::info,  msg); }
    void warn (const std::string &msg) { log(LogLevel::warn,  msg); }
    void error(const std::string &msg) { log(LogLevel::error, msg); }

    /* Returns the process-lifetime default logger (writes to stdout/stderr). */
    static Logger &defaultLogger();
};

/* Default implementation: info→stdout, warn/error→stderr, debug suppressed. */
class StdioLogger : public Logger {
public:
    void log(LogLevel level, const std::string &msg) override;
};

#endif // LOGGER_HPP
