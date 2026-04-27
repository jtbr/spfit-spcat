#ifndef DPI_CONTEXT_HPP
#define DPI_CONTEXT_HPP

#ifdef __cplusplus
extern "C" {
#endif

struct DpiContext
{
  int nvib;
  int iwhole;
  int isdgn;
  int nqn;
};

#ifdef __cplusplus
} // extern "C"
#endif

#endif // DPI_CONTEXT_HPP
