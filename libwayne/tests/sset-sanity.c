#include "misc.h"
#include "sets.h"
#include "rand48.h"
#include <stdio.h>
#include <assert.h>

#define SETSIZE 5

int main()
{
    int n = SETSIZE, i;
    TSET A, B; // don't need to allocate TSETs or SSETs
    SetStartup();
    TSetReset(A);
    TSetReset(B);

    srand48(time(0) + getpid());

    printf("setA setB\n");
    for(i=0; i<n; i++)
    {
	int newA = lrand48() % SETSIZE, newB = lrand48() % SETSIZE;
	TSetAdd(A, newA);
	TSetAdd(B, newB);
	printf("%3d %3d\n", newA, newB);
    }

    puts("Checking sanity...");
    assert(TSetCardinality(A) <= SETSIZE);
    assert(TSetCardinality(B) <= SETSIZE);

    /* N(A & B ) == N(A) + N(B) - N(A union B) */
    assert(TSetCardinality(TSetIntersect(A, B)) ==
	TSetCardinality(A) + TSetCardinality(B) - TSetCardinality(TSetUnion(A, B))
	);
    puts("The world makes sense!");
}
