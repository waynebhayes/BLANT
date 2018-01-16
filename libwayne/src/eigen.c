#include "eigen.h"

int Eigen(int N, double A[N][N], double eigval[N], double eigvec[N][N])
{
    const int LWORK = 3*N-1;
    double WORK[LWORK];
    int i, j, INFO;
    char jobz, uplo;
    extern void dsyev_();

    jobz='V';
    uplo='L';   /* upper part in C, but F sees it transposed */
    
    for(i=0; i<N; i++)
	for(j=i; j<N; j++)
	    eigvec[i][j] = A[i][j];

    dsyev_(&jobz, &uplo, &N, eigvec, &N, eigval, WORK, &LWORK, &INFO);

    if(INFO < 0)
	Warning("Eigen: argument %d illegal input", -INFO);
    else if(INFO > 0)
	Warning("Eigen: failed to converge, %d off-diag elem", INFO);

    return INFO;
}
