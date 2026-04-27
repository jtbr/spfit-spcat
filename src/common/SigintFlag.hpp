/*   Copyright (C) 1989, California Institute of Technology */
/*   All rights reserved.  U. S. Government Sponsorship under */
/*   NASA Contract NAS7-918 is acknowledged. */

#ifndef SIGINTFLAG_HPP
#define SIGINTFLAG_HPP

#include <csignal>

/* RAII replacement for rqexit/brkqr: installs SIGINT handler on construction,
   exposes isTriggered(), restores prior handler on destruction. */
class SigintFlag {
public:
    SigintFlag();
    ~SigintFlag();
    static bool isTriggered();

private:
    static void handler(int);
    static volatile sig_atomic_t s_flag;
    void (*m_prev_handler)(int);
};

/* chrono-based replacement for caldelay(): returns true once per interval.
   Resets the interval when delay_seconds changes. */
bool caldelay(int delay_seconds);

#endif // SIGINTFLAG_HPP
