#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "misc.h"
#include "matvec.h"
#include "f_eval.h"

typedef struct _rk23
{
    int n;
    double TOL,HMAX,HMIN,H,T,*w;
    F_EVAL F;
} RK23;


RK23 *Rk23Alloc(int n, double t, double *y, F_EVAL f, int stiff_flag,
    double dt, double zero);

double Rk23Integrate(RK23 *r, double B);

#define Rk23Free(r) Free(r)
