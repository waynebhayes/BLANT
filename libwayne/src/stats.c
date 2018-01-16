/* Version 0.0
** From "Wayne's Little DSA Library" (DSA == Data Structures and
** Algorithms) Feel free to change, modify, or burn these sources, but if
** you modify them please don't distribute your changes without a clear
** indication that you've done so.  If you think your change is spiffy,
** send it to me and maybe I'll include it in the next release.
**
** Wayne Hayes, wayne@cs.utoronto.ca (preffered), or wayne@cs.toronto.edu
*/

#include "misc.h"
#include "stats.h"
#include "rand48.h"
#include <math.h>

STAT *StatAlloc(int numHistogramBins, double histMin, double histMax,
    Boolean geom)
{
    STAT *s = Calloc(1, sizeof(STAT));
    s->n = 0;
    s->geom = geom;
    s->sum = s->sum2 = s->sum3 = s->geomSum = s->geomSum2 = s->geomSum3 = 0.;
    s->min = 1e30;
    s->max = -1e30;
    if(numHistogramBins)
    {
	s->histogram = Calloc(1, (numHistogramBins+2) * sizeof(int));
	++s->histogram; /* make the base -1 */
	s->numHistBins = numHistogramBins;
	s->histMin = histMin;
	s->histWidth = histMax-histMin;
	s->histCumulative = false;
    }
    return s;
}

STAT *StatReset(STAT *s)
{
    s->n = 0;
    s->sum = s->sum2 = s->sum3 = s->geomSum = s->geomSum2 = s->geomSum3 = 0.;
    if(s->numHistBins)
    {
	memset(s->histogram - 1, 0, (s->numHistBins+2) * sizeof(int));
	s->histCumulative = false;
    }
    return s;
}

void StatFree(STAT *s)
{
    if(s->numHistBins)
	free(s->histogram-1);   /* since we bumped it up before */
    free(s);
}


static void ToggleHistType(STAT *s)
{
    int i;
    if(!s->histCumulative)	/* make it cumulative */
    {
	s->histCumulative = true;
	for(i=0; i <= s->numHistBins; i++)
	    s->histogram[i] += s->histogram[i-1];
    }
    else
    {
	s->histCumulative = false;
	for(i=s->numHistBins; i >= 0; i--)
	    s->histogram[i] -= s->histogram[i-1];
    }
}


void StatAddSample(STAT *s, double sample)
{
    s->n++;
    s->sum += sample;
    s->sum2 += sample * sample;
    s->sum3 += sample * sample * sample;
    if(sample > s->max)
	s->max = sample;
    if(sample < s->min)
	s->min = sample;
    if(s->geom)
    {
	if(sample > 0)
	{
	    double ls = log(sample);
	    s->geomSum += ls;
	    s->geomSum2 += ls*ls;
	    s->geomSum3 += ls*ls*ls;
	}
	else
	    Warning("StatAddSample(geom): sample %g <= 0", sample);
    }
    if(s->numHistBins)
    {
	int histBin = s->numHistBins *
	    (sample - s->histMin)/s->histWidth;
	if(s->histCumulative)
	    ToggleHistType(s);
	if(histBin < 0)
	    ++s->histogram[-1];
	else if(histBin >= s->numHistBins)
	    ++s->histogram[s->numHistBins];
	else
	    ++s->histogram[histBin];
    }
}

void StatDelSample(STAT *s, double sample)
{
    s->n--;
    s->sum -= sample;
    s->sum2 -= sample*sample;
    s->sum3 -= sample * sample * sample;
    if(sample <= s->min)
    {
	if(sample < s->min)
	    Warning("StatDelSample: deleted sample is less than the observed minimum!");
	else
	    Warning("StatDelSample: can't update minimum when it's the deleted sample");
    }
    if(sample >= s->max)
    {
	if(sample > s->max)
	    Warning("StatDelSample: deleted sample is greater than the observed maximum!");
	else
	    Warning("StatDelSample: can't update maximum when it's the deleted sample");
    }
    if(s->geom)
    {
	if(sample > 0)
	{
	    double ls = log(sample);
	    s->geomSum -= ls;
	    s->geomSum -= ls*ls;
	    s->geomSum -= ls*ls*ls;
	}
	else
	    Warning("StatDelSample(geom): sample %g <= 0", sample);
    }
    if(s->numHistBins)
    {
	int histBin = s->numHistBins *
	    (sample - s->histMin)/s->histWidth;
	if(s->histCumulative)
	    ToggleHistType(s);
	if(histBin < 0)
	    --s->histogram[-1];
	else if(histBin >= s->numHistBins)
	    --s->histogram[s->numHistBins];
	else
	    --s->histogram[histBin];
    }
}

int *StatHistogram(STAT*s)
{
    if(s->histCumulative)
	ToggleHistType(s);
    return s->histogram;
}

int *StatCumulativeHistogram(STAT*s)
{
    if(!s->histCumulative)
	ToggleHistType(s);
    return s->histogram;
}


double StatMean(STAT*s)
{
    return s->sum / s->n;
}

double StatGeomMean(STAT*s)
{
    return exp(s->geomSum / s->n);
}

double StatVariance(STAT*s)
{
    /* unbiased estimator, Law & Kelton eqn 4.4 */
    return (s->sum2 - s->sum*s->sum / s->n) / (s->n - 1);
}

double StatGeomVariance(STAT*s)
{
    return exp((s->geomSum2 - s->geomSum*s->geomSum / s->n) / (s->n - 1));
}

double StatStdDev(STAT*s)
{
    return sqrt(StatVariance(s));
}

double StatGeomStdDev(STAT*s)
{
    return exp(sqrt(log(StatGeomVariance(s))));
}

double StatSkew(STAT*s)
{
    /* this is probably not unbiased; Law & Kelton p.290.  I've substituted
    ** Xbar for mu (tsk, tsk), and expanded the E[(x-mu)^3] so I can use
    ** running sums.
    */
    double mu = s->sum / s->n;
    double sd = StatStdDev(s);
    return ((s->sum3 + 3*mu*s->sum2)/s->n + 4*mu*mu*mu) / (sd*sd*sd);
}

double StatConfInterval(STAT *s, double confidence)
{
    return StatTDist((1-confidence)/2, s->n - 1) * sqrt(StatVariance(s) / s->n);
}


double StatRV_Normal(void)
{
    static int which;
    static double next;
    double fac, rsq, v1, v2;

    if(!which)
    {
	do {
	    v1 = 2*drand48()-1;
	    v2 = 2*drand48()-1;
	    rsq = SQR(v1)+SQR(v2);
	} while(rsq >= 1 || rsq == 0);
	fac=sqrt(-2*log(rsq)/rsq);
	next = v1*fac;
	which = 1;
	return v2*fac;
    }
    else
    {
	which = 0;
	return next;
    }
}


/* Taken from MacDougall, "Simulating Computer * Systems", MIT Press, 1987,
 * p. 276 (or therabouts.)
 */
static double NormalDist (double quantile)
{
    double    q, z1, n, d;

    q = quantile > 0.5 ? (1 - quantile) : quantile;
    z1 = sqrt (-2.0 * log (q));
    n = (0.010328 * z1 + 0.802853) * z1 + 2.515517;
    d = ((0.001308 * z1 + 0.189269) * z1 + 1.43278) * z1 + 1.0;
    z1 -= n / d;
    return (quantile > 0.5 ? -z1 : z1);
}


/*
 * TDist (double quantile, long freedom) computes the given upper quantile
 * of the student t distribution with "freedom" degrees of freedom.  This
 * is the x for which the area under the curve from x to +ve infinity is
 * equal to quantile.  Taken from MacDougall, "Simulating Computer
 * Systems", MIT Press, 1987, p. 276.
 */

double StatTDist (double quantile, long freedom)
{
    long    i;
    double    z1, z2, h[4], x;

    z1 = fabs (NormalDist (quantile));
    z2 = z1 * z1;

    h[0] = 0.25 * z1 * (z2 + 1.0);
    h[1] = 0.010416667 * z1 * ((5.0 * z2 + 16.0) * z2 + 3.0);
    h[2] = 0.002604167 * z1 * (((3.0 * z2 + 19.0) * z2 + 17.0) * z2 - 15.0);
    h[3] = z1 * ((((79.0 * z2 + 776.0) * z2 + 1482.0) * z2 - 1920.0)
	 * z2 - 945.0);
    h[3] *= 0.000010851;

    x = 0.0;
    for (i = 3; i >= 0; i--)
    x = (x + h[i]) / (double) freedom;
    z1 += x;
    return (quantile > 0.5 ? -z1 : z1);
}

#if TEST
main()
{
    int v;
    double gamma;
    puts("Enter df, gamma pairs until you're happy (see Law&Kelton, Appendix)");
    while(scanf("%d %lf", &v, &gamma) == 2)
	printf("%g\n", StatTDist(gamma, v));
}
#endif
