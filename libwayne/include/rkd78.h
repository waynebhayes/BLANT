#ifndef _RKD78_H
#define _RKD78_H

#include "f_eval.h"

typedef struct _rkd78
{
    double T, *Y, TOL, *C, *W;
    int N, IND, NW;
    F_EVAL F;
} RKD78;

/*
** RKD78 is a 7/8 pair Runge-Kutta method using Defect-based error
** control, written by Wayne Enright (+ others?)
** n = number of equations.
** time = initial time.
** y = pointer to your array of n state variables.
** f = function of yours that computes ydot
** eps = allowable relative error.
**
** stiff_flag must be 0, and zero must 0.0.
*/ 
RKD78 *Rkd78Alloc(int n, double time, double *y, F_EVAL f,
    int stiff_flag, double eps, double zero);

/* Do the actual integration.  Return the actual TOUT */
double Rkd78Integrate(RKD78 *r, double TOUT);

void Rkd78Free(RKD78 *r);

#endif  /* _RKD78_H */
