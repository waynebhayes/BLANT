#ifndef _LONGLONG_H
#define _LONGLONG_H

/* These routines should handle 64 bit signed long longs.
*/
#include <stdio.h>

char *lltoa(char *foo, long long n);	/* char* can be NULL */
#ifndef __sun__
//long long atoll(char *p);
#endif
long long getll(FILE *fp);

#define LONG_LONG_BITS 64
#define MAX_ASCII_LONG_LONG 21
#define MAX_LLTOA_BUFS 100

#endif  /* _LONGLONG_H */
