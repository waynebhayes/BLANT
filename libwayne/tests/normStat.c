#include "misc.h"
#include "stats.h"

int main(int argc, char *argv[])
{
    int i, N = atoi(argv[1]), M = atoi(argv[2]);

    srand48(time(0)+getpid());
    for(i=0; i<N; i++)
    {
	int j;
	double V = 0;
	for(j=0; j<M; j++)
	{
	    double Z = StatRV_Normal(), U = SQR(Z);
	    V += U;
	}
	printf("%g\n", V);
    }

    return 0;
}
