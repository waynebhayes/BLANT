#include "misc.h"
#include "stats.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
    PEARSON *p;
    FILE *fp = NULL;
    int nextArg=1;

    p = PearsonAlloc();

    assert(nextArg <= argc);
    if(nextArg == argc) fp = stdin;
    else fp = Fopen(argv[nextArg], "r");

    double x,y;
    while(fscanf(fp, "%lf %lf", &x, &y) == 2)
	PearsonAddSample(p, x,y);

    puts(PearsonPrint(p));

    PearsonFree(p);

    return 0;
}
