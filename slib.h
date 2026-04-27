#ifndef SLIB_H
#define SLIB_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/************** SLIB interface ***************************************** */
/* Implementations live in ulib.c (moved from slibgcc.c). */
int  chtime(char str[], const int n);
int  fgetstr(/*@out@*/char buffer[], const int n, FILE *stream);

#ifdef __cplusplus
}
#endif

#endif /* SLIB_H */
