#ifndef _ADAPTIVE_VERLET_H
#define _ADAPTIVE_VERLET_H

/*
** Note that the function only takes the positions, not the velocities,
** as input; and it returns the accelerations, only.
*/
typedef void (*FF_EVAL)(int N, double T, double *R, double *Rdotdot);

void init_adaptive_verlet(int n, double t0, double dt, double *r, double *v, FF_EVAL f);

/*
** return actual tout.
*/
double integrate_adaptive_verlet(double tout);

#endif  /* _ADAPTIVE_VERLET_H */
