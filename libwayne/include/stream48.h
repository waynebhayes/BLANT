/*  math.c  --  should be in global math library  */
#include <math.h>
#include <unistd.h>

/***********************************************************************
**  rand48 multiple stream code. You must set _rand48streams BEFORE
** you call stream48() if you want more than a single stream. Thereafer
** just use stream48() to choose which stream you want to affect before
** calling the function you want. For example, if you want 5 streams, you
** set _rand48streams to 5 somewhere before calling any of the rand48
** functions. Then you could do stuff like:
**
**  stream48 (0);       // the default stream ;
**  srand48 (4345);     // seed stream 0 with 4345
**  stream48 (1);
**  srand48 (6);        // give stream 1 a different seed ...
**
**  etc.
** Stream() returns the previous stream. (the one you just replaced).
**********************************************************************/

extern int _rand48streams;	/* default 1 */

/* Must be setbefore ANY of the routines below are called. */
#define Stream48Init(n) _rand48streams = (n)
int Stream48Which(void);  // return current stream number.

int Stream48(int n);	/* choose the current stream */

/* Some helpers. */
long Stream48Randomize(void);	 /* randomize current stream.  Return seed. */
long Stream48RandInt(long minimum, long maximum); /* min, max inclusive */
double Stream48Rand(void);	/* drand48 with min of 1e-10 */

