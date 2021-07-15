#include "misc.h"
#include "stats.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
    FILE *fp = NULL;
    double x,y;
    int nextArg=1;

    COVAR *c = CovarAlloc();

    assert(nextArg <= argc);
    if(nextArg == argc) fp = stdin;
    else fp = Fopen(argv[nextArg], "r");

    while(fscanf(fp, "%lf%lf", &x,&y) == 2)
	CovarAddSample(c,x,y);

    printf("# %6d covar %.6g\n", c->n, Covariance(c));

    CovarFree(c);

    return 0;
}
