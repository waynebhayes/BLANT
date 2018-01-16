#ifndef _COMPRESSEDINT_H
#define _COMPRESSEDINT_H
#include <stdio.h>

long CompressedIntRead(FILE *fp);
long CompressedIntWrite(long value, FILE *fp);
#endif  /* _COMPRESSEDINT_H */
