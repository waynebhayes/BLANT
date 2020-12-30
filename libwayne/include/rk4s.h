#ifndef _RK4S_H
#define _RK4S_H

#include "misc.h"
#include "f_eval.h"

typedef struct _rk4s
{
    int N;
    Boolean Hamiltonian;
    double DT, *userY, *internalY, T;
    F_EVAL F;
} RK4S;


/*
 * RK4S is a symplectic 5-stage, 4th order RK formula from
 * _Numerical_Hamiltonian_Problems_ by Sanz-Serna & Calvo (1994).
 */
RK4S *Rk4sAlloc(int n, double t, double *y, F_EVAL f, int stiff_flag,
    double dt, double zero);

/* Do the actual integration.  Return the actual TOUT */
double Rk4sIntegrate(RK4S *l, double TOUT);

void Rk4sFree(RK4S *l);

#endif  /* _RK4S_H */
