#ifndef _IQRKD78_H
#define _IQRKD78_H

#include "f_eval.h"

typedef struct _iqrkd78
{
    double T, *Y, TOL;
    long double QT, *QY, QTOL, *C, *W;
    int N, IND, NW;
    F_EVAL F;
} IQRKD78;

/*
** RKD78 is a 7/8 pair Runge-Kutta method using Defect-based error
** control, written by Wayne Enright (+ others?)  QRKD is a
** quad-precision vesion of it.  Fortran's QUAD==gcc's long double.
** Finally, IQRKD is a double-prec interface to the quad routines,
** for allowing more accurate integrations even though it looks like
** a double-precision routine from the outside.
** n = number of equations.
** time = initial time.
** y = pointer to your array of n state variables.
** f = function of yours that computes ydot
** eps = allowable relative error.
**
** stiff_flag must be 0, and zero must 0.0.
*/
IQRKD78 *Iqrkd78Alloc(int n, double t0, double *y, F_EVAL f,
    int stiff_flag, double eps, double zero);

/* Do the actual integration.  Return the actual TOUT */
double Iqrkd78Integrate(IQRKD78 *r, double TOUT);

void Iqrkd78Free(IQRKD78 *r);

#endif  /* _IQRKD78_H */
