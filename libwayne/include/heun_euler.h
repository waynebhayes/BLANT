#ifndef _HEUN_EULER_H
#define _HEUN_EULER_H

/*
** Note that the function only takes the positions, not the velocities,
** as input; and it returns the accelerations, only.
*/
typedef void (*FF_EVAL)(int N, double T, double *R, double *Rdotdot);

// set tol to 0 or negative for constant timestep.
void init_heun_euler(int n, double t0, double dt, double *r, double *v, FF_EVAL f, double tol);

/*
** return actual tout.
*/
double integrate_heun_euler(double tout);

#endif  /* _HEUN_EULER_H */
