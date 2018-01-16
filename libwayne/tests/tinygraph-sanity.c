#include <stdio.h>
#include "misc.h"
#include "tinygraph.h"

int main(int argc, char *argv[])
{
    int BFSsize, i, j, n;
    TINY_GRAPH *G, *Gbar, *GG;
    for(n=100000;n>=0;n--)
    {
	G->n = 6 + lrand48() % 3; // ensure G->n is at least 2
	G = TinyGraphAlloc(G->n);
	for(int k=2*G->n; k>=0; k--)
	{
	    j = i = lrand48() % G->n;
	    while(j==i) j = lrand48() % G->n;
	    TinyGraphConnect(G,i,j);
	}
	Gbar = TinyGraphComplement(NULL, G);
	GG = TinyGraphComplement(NULL, Gbar);
	for(i=0; i<G->n; i++)
	    assert(G->A[i] == GG->A[i]);
	int perm[G->n];
	if(TinyGraphsIsomorphic(perm, G, Gbar))
	{
	    TinyGraphPrintAdjMatrix(stdout, G);
	    printf("G(%d) is isomorphic to Gbar, perm = [", G->n);
	    for(i=0; i<G->n; i++) printf(" %d", perm[i]);
	    puts(" ]");
	}
	TinyGraphFree(G);
	TinyGraphFree(GG);
	TinyGraphFree(Gbar);
    }

    while(false && !feof(stdin))
    {
	int root, distance, nodeArray[MAX_TSET], distArray[MAX_TSET];
	printf("enter BSF start node, distance: "); fflush(stdout);
	if(scanf("%d%d", &root, &distance) != 2)
	    break;
	BFSsize = TinyGraphBFS(G, root, distance, nodeArray, distArray);
	for(i=0; i<BFSsize; i++)
	    printf("%d(%d) ", nodeArray[i], distArray[nodeArray[i]]);
	printf("\n");
    }

    return 0;
}
