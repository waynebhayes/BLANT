#ifndef _LEAPFROG_H
#define _LEAPFROG_H

/*
** Note that the function only takes the positions, not the velocities,
** as input; and it returns the accelerations, only.
*/
typedef void (*FF_EVAL)(int N, double T, double *R, double *Rdotdot);

void init_leapfrog(int n, double t0, double dt, double *r, double *v, FF_EVAL f);

/*
** return actual tout.
*/
double integrate_leapfrog(double tout);

#endif  /* _LEAPFROG_H */
