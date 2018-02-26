#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "misc.h"
#include "stats.h"

char USAGE[] = "USAGE: %s batchSize precision confidence\n"
    "read numbers on the standard input (no files!), in batches of batchSize\n"
    "keep reading until we know the mean within precision (smaller is better),\n"
    "with confidence 'confidence' (closer to 1 is better)\n";

int main(int argc, char *argv[])
{
    STAT *batch, *batchMeans;
    int batchSize = atoi(argv[1]), numBins=0;
    double sample, precision=atof(argv[2]), confidence=atof(argv[3]), histMin=0, histMax=0;
    Boolean geom = false;

    batch = StatAlloc(numBins, histMin, histMax, geom);
    batchMeans = StatAlloc(numBins, histMin, histMax, geom);

    while(scanf("%lf", &sample) == 1 && (StatSampleSize(batchMeans)<3 || fabs(StatConfInterval(batchMeans, confidence)) > precision))
    {
	StatAddSample(batch, sample);
	if(StatSampleSize(batch) == batchSize)
	{
	    StatAddSample(batchMeans, StatMean(batch));
	    StatReset(batch);
	}
    }

    StatFree(batch);
    printf("# %d mean %.16g min %.16g max %.16g stdDev %.16g var %.16g skew %.16g\n",
	StatSampleSize(batchMeans), StatMean(batchMeans), StatMin(batchMeans), StatMax(batchMeans),
	    StatStdDev(batchMeans), StatVariance(batchMeans), StatSkew(batchMeans));
    StatFree(batchMeans);

    return 0;
}
