#ifndef _MISC_H
#define _MISC_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h> // because I use the time() function so much

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define ABS(a) ((a)>=0?(a):-(a))
#define SQR(x) ((x)*(x))
#define SIGN(x) ((x)==0?0:(x)<0?-1:1)
#define SIGN2(x,y) ((x)*SIGN(y))
#define IMPLIES(a,b) (!(a)||(b))
static __inline__ double sqr(double x) { return x*x; }
static __inline__ double cube(double x) { return x*x*x; }
#define MACH_EPS 1e-16  /* something that really should be in math.h */


/* FA accesses a 1D FORTRAN array using FORTRAN semantics, ie, from 1
 * If you compiled your FORTRAN code using -r8, then integers take up
 * 64 bits, although half of them are not used. (Repeat after me:
 * "I love FORTRAN, honest I do!")
 */

#ifdef sgi
#define R8 1    /* 1 for no, 2 for yes if you compiled f77 -r8 */
#else
#define R8 1
#endif

#define IFA(name, index) (*((name) + R8*((index)-1)))   /* integer array */
#define RFA(name, index) (*((name) + (index)-1))        /* Real array */

#define elsif else if
#define until(x) while(!(x))

#if 0
extern int fprintf(const FILE *, const char *, ...);
extern int printf(const char *, ...);
extern char *puts(char[]), fputs(char[], FILE*);
extern int fflush(FILE*);
extern int fscanf(const FILE *, const char *, ...);
extern int scanf(const char *, ...);
#endif

typedef unsigned char Boolean;
/* enum _boolEnum { false=0, true=1, maybe=2 };*/
#ifndef false
#define false (Boolean)0
#define true  (Boolean)1
#endif
#define maybe (Boolean)2

typedef union _voidInt {
    long long ll;
    long int l, l2;
    unsigned long long ull;
    void *v, *v2;
    char *s, *s2;
    char c[sizeof(long long)];
    int i, i2;
    unsigned int ui, ui2;
    float f, f2;
    double d;
} foint;

/* The comparison function type: used by heaps, binary trees and sorts.
*/
typedef int (*pCmpFcn)(foint, foint);   

/* Copy a foint.  In all instances, you are expected to know what the
** foint actually is, and return a copy of it.  If a FointCopy function
** pointer is ever NULL, the code will do a shallow copy.
*/
typedef foint (*pFointCopyFcn)(foint);   
typedef void (*pFointFreeFcn)(foint);

/* this is is the general error return value for most abstract data types */
extern const foint ABSTRACT_ERROR;

/*
** The following three functions print out warnings, apologies (ie,
** fatal limitations in the program that could be worked around with
** more effort on the part of the library creator), and fatal errors.
** The messages are always sent to stderr.  In addition, it is assumed
** that these messages should also *always* be printed to the terminal,
** so if neither stdout nor stderr are a terminal, then we open /dev/tty
** and put it there anyway.  None of these touch stdout.
*/
extern void Warning(const char *fmt, ...);
extern void Apology(const char *fmt, ...);
extern void Abort(const char *fmt, ...);  /* generates an assertion failure */
extern void Fatal(const char *fmt, ...);  /* Fatal user error */

// Try to use mmap? As long as mmap is supported, even if it doesn't work blant will still work
// because it'll revert to simply loading the entire binary mapping if the mmap fails.
#define MMAP 1
extern void *Mmap(void *p, size_t n, int fd);
extern void *Malloc(size_t);
extern void *Calloc(size_t, size_t);
extern void *Realloc(void *ptr, size_t newSize);
void *Memdup(void *v, size_t n);
#define Free(s) free((void*)(s))


/* return current user time used in seconds.  Returns -1 if error. */
double uTime(void);

char *Int2BitString(char word[33], unsigned i);

/*
** Be careful!  These may not work as you think if i is not unsigned.
** Rotate s bits of i by n positions and mask out anything above bit s.
*/
#if __sun__ /* gcc 2.6.3 bug */
unsigned long long _NotZero;
#define RotLeft(type,i,s,n) ((_NotZero = ~(type)0), ((((i) << (n)) | ((i) >> (s-(n)))) & (_NotZero >> (8*sizeof(i)-(s)))))
#define RotRight(type,i,s,n) ((_NotZero = ~(type)0), ((((i) >> (n)) | ((i) << (s-(n)))) & (_NotZero >> (8*sizeof(i)-(s)))))
#else
#define RotLeft(type,i,s,n) ((((i) << (n)) | ((i) >> (s-(n)))) & (~(type)0 >> (8*sizeof(i)-(s))))
#define RotRight(type,i,s,n) ((((i) >> (n)) | ((i) << (s-(n)))) & (~(type)0 >> (8*sizeof(i)-(s))))
#endif

void PrintArray(FILE *fp, int n, int *array);

long long IIntPow(int, int);
double IntPow(double, int);
int Log2(int);	/* integer part of log to the base 2 */

/* Greatest Common Divisor */
int gcd(int a, int b);
int IsPrime(long long n);
int PrimeFactors(int n, int count[]);

#endif  /* _MISC_H */
