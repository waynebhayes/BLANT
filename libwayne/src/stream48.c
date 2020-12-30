/*  math.c  --  should be in global math library  */
#include "misc.h"
#include "rand48.h"
#include "stream48.h"

long Stream48Randomize (void)
{
    static long seed;
    seed += time((long *) 0)+getpid();
    srand48 (seed);
    return seed;
}


long Stream48RandInt (long minimum, long maximum)
{
    return (long) (drand48 () * (maximum - minimum + 1) + minimum);
}


double Stream48Rand (void)
{
    double  r;

    r = drand48 ();
    return (r == 0.0) ? 1e-10 : r;
}


static unsigned short *equ48 (unsigned short x[3], unsigned short y[3])
{
    x[0] = y[0], x[1] = y[1], x[2] = y[2];
    return x;
}

int _rand48streams = 1;    /* number of streams you want. Must be set
		** before ANY of the routines are called.
		*/

int Stream48 (int n)
{
    static unsigned short (*stream)[3];
    static long whichStream, old;
    static char notReset = 1;
    unsigned short *seed48(unsigned short x[3]);

    if(notReset)
    {
	long i;

	stream = (unsigned short (*)[3])
	    calloc (_rand48streams, sizeof(*stream));

	if(stream == 0)
	{
	    fprintf(stderr, "out of memory\n");
	    exit(1);
	}

	equ48(stream[0], seed48(stream[0]));
	seed48(stream[0]);

	for(i=1; i< _rand48streams; i++)
	    equ48(stream[i], stream[0]);

	notReset = 0;
    }

    if(n < 0 || n >= _rand48streams)
    {
	fprintf(stderr, "rand48stream(): no such stream %d\n", n);
	exit(1);
    }

    if(n != whichStream)
    {
	equ48(stream[whichStream], seed48(stream[n]));
	old = whichStream;
	whichStream = n;
    }

    return old;
}
