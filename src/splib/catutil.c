#include "splib/catutil.h"
int readqn(qnstr, iqn, n)
const char *qnstr;
short *iqn;
const int n;
{ /* read n quanta from string to iqn */
  const char *pqn;
  int i, ich, ic, ival;

  pqn = qnstr;
  i = 0; ival = 0;
  do {
    iqn[i] = 0; ival = 0;
    if (*pqn == 0)
      break;
    ich = (*pqn & 0xff);
    ++pqn;
    if (*pqn == 0)
      break;
    ic = (*pqn & 0xff);
    ++pqn;
    if (ic == ' ' ) 
      continue;
    ival = ic - '0';
    if (ival < 0 || ival > 9)  
      break;
    if (ich != ' ') {
      if (ich == '-') {
        ival = -ival;
      } else if (ich >= '0' && ich <= '9') {
        ival += (ich - '0') * 10;
      } else if (ich >= 'a' && ich <= 'z') {
        ival = -ival - 10 * (ich - ('a' - 1));
      } else if (ich >= 'A' && ich <= 'Z') {
        ival += 10 * (ich - ('A' - 10));
      } else {
        ival = 36;
      }
    }
    iqn[i] = (short) ival;
  } while (++i < n);
  for(ic = i; ic < n; ++ic)
    iqn[ic] = 0;
  if (n == 0 && i == 1) {
    ich = (*pqn & 0xff) - '0';
    if (ich >= 0 && ich <= 9 && ival >= 0) {
      i = 0;
      ival = 10 * ival + ich;
      iqn[0] = (short) ival;
    }
  }
  return (n - i);
}                               /* readqn */

void gupfmt(int igup, char *sgup)
{
  static int czero = (int) '0';
  int i1, i2;
  sgup[0] = ' '; sgup[1] = ' ';
  if (igup <= 9) {
    sgup[2] = (char) (igup + czero);
    return;
  }
  i1 = igup / 10; igup -= i1 * 10; 
  sgup[2] = (char) (igup + czero);
  i2 = i1 / 10; i1 -= i2 * 10;
  sgup[1] = (char) (i1 + czero);
  if (i2 == 0) return;
  if (i2 <= 9) {
    sgup[0] = (char)(i2 + czero);
  } else if(i2 < 36) {
    sgup[0] = (char) (i2 + ((int) 'A' - 10));
  } else {
    sgup[0] = '*';
  }
}
