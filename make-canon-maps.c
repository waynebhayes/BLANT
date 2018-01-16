#include <stdio.h>
#include "misc.h"
#include "tinygraph.h"
#include "blant.h"

// We use TinySets since they can be stored in a char (8 bits), and thus the adjacency matrix for k<=8 fits into 8 chars.
static TINY_GRAPH *G, *_canonicalGraph[MAX_CANONICALS];	// G is reused for each value of numEdges
static int k, _numCanonicals, _canonicalSig[MAX_CANONICALS];


/* The array perm[] will hold the permutation from the canonical to the non-canonical.
** More formally, if perm[i]=j, it means that you should map node i from the canonical
** to node j of the non-canonical in order to prove the isomorphism. Then, when faye is doing
** samples, it means that if we print out the sampled nodes in the order listed in the perm
** array, then columns representing the same canonical present a perfect local alignment.
*/
static int perm[MAX_TSET];

/*
** Given an integer, build the graph into the TINY_GRAPH *G, which has already been allocated.
** Handles either upper or lower triangle representation depending upon compile-time option below.
*/
static void BuildGraph(int Gint)
{
    int i, j, bitPos=0;
    int Gint2 = Gint;  // Gint2 has bits nuked as they're used, so when it's zero we can stop.
    TinyGraphEdgesAllDelete(G);
#if LOWER_TRIANGLE
    for(i=k-1;i>0;i--)
    {
	for(j=i-1;j>=0;j--)
#else	// UPPER_TRIANGLE
    for(i=k-2;i>=0;i--)
    {
	for(j=k-1;j>i;j--)
#endif
	{
	    if(!Gint2) break;
	    int bit = (1 << bitPos);
	    if(Gint & bit)
		TinyGraphConnect(G,i,j);
	    Gint2 &= ~bit;
	    bitPos++;
	}
	if(!Gint2) break;
    }
}


/*
** Given an integer, build the graph and see if it's isomorphic to any existing canonicals.
** If not, then it becomes the next new canonical.
*/
void CheckGraph(int Gint)
{
    BuildGraph(Gint);
    int i, j;
    for(i=0; i<_numCanonicals; i++)
    {
#if PERMS_CAN2NON // the permutation from canonical to non-canonical, which is preferred so columns map to alignments
	if(TinyGraphsIsomorphic(perm, _canonicalGraph[i], G))
#else  // permutation from non-canonical to canonical
	if(TinyGraphsIsomorphic(perm, G, _canonicalGraph[i]))
#endif
	{
	    printf("%d\t%d\t", Gint, _canonicalSig[i]); // We found its canonical, print both
	    for(j=0; j<k; j++) printf("%d", perm[j]);
	    puts("");
	    break;
	}
    }
    if(i == _numCanonicals) // it is a new canonical
    {
	printf("%d\t%d\t", Gint, Gint); // It's a new canonical, so print its own value twice.
	    for(j=0; j<k; j++) printf("%d", j);	// and its permutation is the identity.
	int nodeArray[k], distArray[k];
	if(TinyGraphBFS(G, 0, k, nodeArray, distArray) == k)
	    puts(" 1"); // connected, thus graphlet
	else
	    puts(" 0"); // disconnected, thus not graphlet
	_canonicalSig[_numCanonicals] = Gint;
	TinyGraphCopy(_canonicalGraph[_numCanonicals], G);
	++_numCanonicals;
    }
}

int main(int argc, char *argv[])
{
    int i;
    if(argc != 4) {
	fprintf(stderr, "USAGE: %s k numParallel coreID\n", argv[0]);
	exit(1);
    }
    k = atoi(argv[1]);
    assert(k > 2 && k <= 8);
    int numParallel=atoi(argv[2]), coreID=atoi(argv[3]);
    G = TinyGraphAlloc(k);
    for(i=0; i<MAX_CANONICALS; i++)
	_canonicalGraph[i] = TinyGraphAlloc(k);

    int maxGint = (1<<(k*(k-1)/2)), Gint;
    int numPerCore = maxGint / numParallel;
    int Gmin = coreID*numPerCore;
    int Gmax = (coreID+1)*numPerCore;
    for(Gint=Gmin; Gint < Gmax; Gint++)
	CheckGraph(Gint);

    return 0;
}
