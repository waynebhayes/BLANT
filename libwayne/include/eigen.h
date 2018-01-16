/*
** Eigen: compute eigenvalues & eigenvectors of SYMMETRIC input matrix A.
** Since the matrix is symmetric, it is sufficient if the upper right of
** A is filled.
** Internally we call LAPACK's dsyev.  Return value is dsyev's INFO:
**
** = 0:  successful exit
** < 0:  if INFO = -i, the i-th argument had an illegal value
** > 0:  if INFO = i, the algorithm failed to converge; i off-diagonal
** elements of an intermediate tridiagonal form did not converge to zero.
**
** row (not column) i of eigvec is the eigenvector corresponding to eigval[i].
*/

#define EigenLogMaxN 8
#define EigenMaxN (1<<EigenLogMaxN)

int Eigen(int N, double A[N][N], double eigval[N], double eigvec[N][N]);
