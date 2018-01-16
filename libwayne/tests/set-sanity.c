#include "misc.h"
#include "sets.h"
#include "rand48.h"
#include <stdio.h>
#include <assert.h>

#define SETSIZE 20

int main()
{
    int n = SETSIZE, i;
    SET *A = SetAlloc(SETSIZE);
    SET *B = SetAlloc(SETSIZE);
    SET *tmp = SetAlloc(SETSIZE);

    srand48(time(0) + getpid());

    printf("setA setB\n");
    for(i=0; i<n; i++)
    {
	int newA = lrand48() % SETSIZE, newB = lrand48() % SETSIZE;
	SetAdd(A, newA);
	SetAdd(B, newB);
	printf("%3d %3d\n", newA, newB);
    }

    /* add a few more, just for kicks */
    SetAddList(A, 7, 4, 2, 9, -1);

    puts("Checking sanity...");
    assert(SetCardinality(A) <= SETSIZE);
    assert(SetCardinality(B) <= SETSIZE);

    /* N(A & B ) == N(A) + N(B) - N(A union B) */
    assert(  SetCardinality(SetIntersect(tmp, A, B)) ==
	SetCardinality(A) + SetCardinality(B) - SetCardinality(SetUnion(tmp, A, B))
	);
    puts("The world makes sense!");
}
