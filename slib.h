#ifndef FILE
#include <stdio.h>
#endif
/************** SLIB interface ***************************************** */
/* Implementations live in ulib.c (moved from slibgcc.c). */
int  chtime(char str[], const int n);
int  fgetstr(/*@out@*/char buffer[], const int n, FILE *stream);
