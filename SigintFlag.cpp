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

bool SigintFlag::triggered() const { return s_flag != 0; }

void SigintFlag::handler(int) { s_flag = 1; }
