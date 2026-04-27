#include <chrono>
#include "SigintFlag.hpp"

volatile sig_atomic_t SigintFlag::s_flag = 0;

SigintFlag::SigintFlag()
{
    s_flag = 0;
    m_prev_handler = signal(SIGINT, handler);
}

SigintFlag::~SigintFlag()
{
    signal(SIGINT, m_prev_handler);
}

bool SigintFlag::isTriggered() { return s_flag != 0; }

void SigintFlag::handler(int) { s_flag = 1; }

bool caldelay(int delay_seconds)
{
    using clock = std::chrono::steady_clock;
    static clock::time_point next_time = clock::now();
    static int last_delay = -1;
    auto now = clock::now();
    if (delay_seconds != last_delay) {
        last_delay = delay_seconds;
        next_time = now + std::chrono::seconds(delay_seconds);
        return delay_seconds <= 0;
    }
    if (now >= next_time) {
        next_time += std::chrono::seconds(delay_seconds);
        return true;
    }
    return false;
}
