#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "misc.h"
#include "matvec.h"
#include "f_eval.h"

typedef struct _rk12
{
    int n;
    double TOL,HMAX,HMIN,H,T,*w;
    F_EVAL F;
} RK12;


RK12 *Rk12Alloc(int n, double t, double *y, F_EVAL f, int stiff_flag,
    double dt, double zero);

double Rk12Integrate(RK12 *r, double B);

#define Rk12Free(r) Free(r)
