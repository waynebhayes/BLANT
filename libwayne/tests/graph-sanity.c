#include <stdio.h>
#include "misc.h"
#include "graph.h"

int main(int argc, char *argv[])
{
    int BFSsize, i, j;
    Boolean sparse=false, supportNames = true;
    GRAPH *G = GraphReadEdgeList(stdin, sparse, supportNames);
    GRAPH *Gbar = GraphComplement(NULL, G);
    GRAPH *GG = GraphComplement(NULL, Gbar);
    printf("Checking sanity of Complement(Complement(G))...");
    assert(GG->n == G->n);
    for(i=0; i<G->n; i++)
    {
	if(!G->sparse) assert(SetEq(GG->A[i], G->A[i]));
	assert(GG->degree[i] == G->degree[i]);
	for(j=0;j<G->n; j++)
	    assert(GG->neighbor[i][j] == GG->neighbor[i][j]);
    }
    puts("passed!");
    printf("Now count connected components via BFS:\n");

    int root, distance, nodeArray[GG->n], distArray[GG->n], CC=0;
    GG = GraphCopy(GG, G);
    Boolean touched[GG->n];
    for(i=0; i<GG->n; i++) touched[i] = false;

    for(root=0; root < GG->n; root++)
    {
	if(!touched[root])
	{
	    BFSsize = GraphBFS(GG, root, GG->n, nodeArray, distArray);
	    printf ("Connected Component %d: BFSsize = %d\n", CC, BFSsize);
	    for(i=0; i<BFSsize; i++)
	    {
		touched[nodeArray[i]] = true;
		//printf("%d(%d) ", nodeArray[i], distArray[nodeArray[i]]);
	    }
	    ++CC;
	}
    }
    printf("Graph has %d connected components\n", CC);

    return 0;
}
