#include "Logger.hpp"
#include <cstdio>

void StdioLogger::log(LogLevel level, const std::string &msg)
{
    if (level == LogLevel::debug) return;
    FILE *out = (level >= LogLevel::warn) ? stderr : stdout;
    fprintf(out, "%s\n", msg.c_str());
}

Logger &Logger::defaultLogger()
{
    static StdioLogger instance;
    return instance;
}
