#ifndef _VARIABLE_LEAPFROG_H
#define _VARIABLE_LEAPFROG_H

#include "leapfrog.h"

/* A user-provided function that is called before every timestep that should
 * return what timestep to use.
 */
typedef double (*TSTEP_EVAL)(int N, double T, double current_DT, double *R, double *V);

/*
** Note that the function only takes the positions, not the velocities,
** as input; and it returns the accelerations, only.
*/

void init_variable_leapfrog(int n, double t0, TSTEP_EVAL dt_eval,
    double *r, double *v, FF_EVAL f);

/*
** return actual tout.
*/
double integrate_variable_leapfrog(double tout);

#endif  /* _VARIABLE_LEAPFROG_H */
