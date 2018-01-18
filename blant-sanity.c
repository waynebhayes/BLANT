// BLANT sanity checker: given a value of k and the name of a "large" graph input file G, read lines from BLANT that
// have been SORTED by the first column, which is the canonical graphlet ID.  Then, without caring which graphlet
// it actually is, simply verify that adjacent lines that have the same value in the first column, contain a list
// of nodes whose column-wise induced subgraph on G are exactly the same graph, node for node and edge for edge.
// (We do *not* fire up an isomorphism checker because BLANT is supposed to have already done that.)
// The output of BLANT that we're reading should be on stdin, so basically you should run us like this:
//
//      ./blant 7 10000 syeast.el | sort -n | ./blant-sanity 7 syeast.el
// We do not need to know anyhting about canonicals or orbits or even whether the representation is upper
// or lower triangular.  Thus no need for "blant.h"

#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"

static int _k;

// Given the big graph G and a set of nodes in V, return the TINY_GRAPH created from the induced subgraph of V on G.
static TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, int *Varray)
{
    unsigned i, j;
    TinyGraphEdgesAllDelete(Gv);
    for(i=0; i < _k; i++) for(j=i+1; j < _k; j++)
        if(GraphAreConnected(G, Varray[i], Varray[j]))
            TinyGraphConnect(Gv, i, j);
    return Gv;
}

// This is the single-core version of blant. It used to be the main() function, which is why it takes (argc,argv[]).
// Now it simply gets called once for each core we use when run with multiple cores.
int main(int argc, char *argv[])
{
    int k=atoi(argv[1]), i, j;
    _k = k;
    FILE *fp = fopen(argv[2], "r");
    GRAPH *G = GraphReadEdgeList(fp, true); // sparse = true
    assert(G->n >= k);
    fclose(fp);

    static char name[BUFSIZ], line[BUFSIZ];
    SET *V = SetAlloc(G->n);
    unsigned Varray[8], a[8];
    while(fgets(line, sizeof(line), stdin))
    {
	int i,j;
	char newName[BUFSIZ];
	switch(k)
	{
	case 3: sscanf(line, "%s %d %d %d ",newName,a,a+1,a+2); break;
	case 4: sscanf(line, "%s %d %d %d %d ",newName,a,a+1,a+2,a+3); break;
	case 5: sscanf(line, "%s %d %d %d %d %d ",newName,a,a+1,a+2,a+3,a+4); break;
	case 6: sscanf(line, "%s %d %d %d %d %d %d ",newName,a,a+1,a+2,a+3,a+4,a+5); break;
	case 7: sscanf(line, "%s %d %d %d %d %d %d %d ",newName,a,a+1,a+2,a+3,a+4,a+5,a+6); break;
	case 8: sscanf(line, "%s %d %d %d %d %d %d %d %d ",newName,a,a+1,a+2,a+3,a+4,a+5,a+6,a+7); break;
	default: Fatal("hmm, unknown value of k %d", k); break;
	}
	if(strcmp(name, newName) != 0)
	{
	    strcpy(name, newName);
	    for(i=0;i<k;i++)Varray[i]=a[i];
	}
	//fputs(line, stdout); fflush(stdout);
	for(i=0; i<k; i++)for(j=i+1;j<k;j++)if(i!=j)
	{
	    //printf("%d-%d, %d-%d\n", a[i], a[j], Varray[i], Varray[j]); fflush(stdout);
	    assert(GraphAreConnected(G, a[i], a[j]) == GraphAreConnected(G, Varray[i], Varray[j]));
	}
    }
    return 0;
}
