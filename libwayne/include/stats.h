#ifndef _STATS_H
#define _STATS_H

/* Cheap statistics taker and computer */

typedef struct _statistic {
    int n, numHistBins;
    Boolean geom, histCumulative;
    double sum, sum2, sum3, geomSum, geomSum2, geomSum3;
    double histMin, histWidth, min, max;
    int *histogram;
} STAT;

/* set numHistogramBins to zero if you don't want a histogram.
** Set geometricStuff to true if you want to keep track of the geometric
** mean.
*/
STAT *StatAlloc(int numHistogramBins, double histMin, double histMax,
    Boolean geometricStuff);
STAT *StatReset(STAT *s);   /* reset the stats of this variable */
void StatFree(STAT *s);     /* de-allocate the variable */
void StatAddSample(STAT *s, double sample);
void StatDelSample(STAT *s, double sample);
double StatMean(STAT*);
#define StatMin(s) ((s)->min)
#define StatMax(s) ((s)->max)
double StatGeomMean(STAT*);
double StatVariance(STAT*);
double StatGeomVariance(STAT*);
double StatStdDev(STAT*);   /* just the square root of the variance */
double StatGeomStdDev(STAT*);
double StatSkew(STAT*);  /* measure of assymetry of the distribution */
double StatGeomSkew(STAT*);
#define StatSampleSize(s) ((s)->n)
int *StatHistogram(STAT*);	/* non-cumulative histogram */
int *StatCumulativeHistogram(STAT*);

/*
** Random number distributions.
*/
double StatRV_Normal(void);	/* return a N(0,1) random variable */

/* only works if histogramming turned on. eg 0.5 for median */
double StatQuantile(STAT*, double quantile);

/*
** Input your desired confidence, it returns the half-interval.  So say
** you want 0.95 (95%) confidence, it will return you a number x meaning
** "we are 95% sure that the interval [StatMean - x, StatMean + x]
** contains the true mean." It assumes the X's are normally distributed;
** if this is horribly false (ie the Skew ** is greater than about 2 or
** 3), and the number of samples is too small, this can be optimistic.
*/
double StatConfInterval(STAT*, double confidence);
double StatTDistP2Z(double quantile, long freedom);

#endif /* _STATS_H */
