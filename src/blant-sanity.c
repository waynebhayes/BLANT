// BLANT sanity checker: given a value of k, n and the name of a "large" graph input file G, read lines from
// BLANT that have been sorted by the first column, which is the canonical graphlet ID.  Then, without caring
// which graphlet it actually is, simply verify that adjacent lines that have the same value in the first column,
// contain of nodes whose column-wise induced subgraph on G are exactly the same graph, node for node and edge
// for edge. a list (We do *not* fire up an isomorphism checker because BLANT is supposed to have already done
// that.) The output of BLANT that we're reading should be on stdin, so basically you should run us like this:
//
//      ./blant 7 10000 syeast.el | sort -n | ./blant-sanity 7 10000 syeast.el

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
#include "graph.h"

static int _k;

#define Gni GraphNodeName2Int

// Get k and the filename of G on the command line. Then read lines on stdin that should be
// sorted lines from the output of BLANT for that value of k. Ensure that adjacent lines that
// are the same canonical graphette are indeed the exact same graph.
int main(int argc, char *argv[])
{
    int k=atoi(argv[1]), n = atoi(argv[2]), lines;
    FILE *fp = fopen(argv[3], "r");
    _k = k;
    Boolean self=false, directed=false, weighted=false;
    GRAPH *G = GraphReadEdgeList(fp, self, directed, weighted);
    assert(G->n >= k);
    fclose(fp);

    lines=0;
    char line[BUFSIZ];
    while(fgets(line, sizeof(line), stdin))
    {
	int i,j, numRead=-1;
	++lines;
	// SET *V = SetAlloc(G->n);
	char newName[BUFSIZ], name[BUFSIZ], Sarray[8][BUFSIZ], sa[8][BUFSIZ];
	switch(k)
	{
	case 3: numRead=sscanf(line, "%s %s %s %s ",newName,sa[0],sa[1],sa[2]); break;
	case 4: numRead=sscanf(line, "%s %s %s %s %s ",newName,sa[0],sa[1],sa[2],sa[3]); break;
	case 5: numRead=sscanf(line, "%s %s %s %s %s %s ",newName,sa[0],sa[1],sa[2],sa[3],sa[4]); break;
	case 6: numRead=sscanf(line, "%s %s %s %s %s %s %s ",newName,sa[0],sa[1],sa[2],sa[3],sa[4],sa[5]); break;
	case 7: numRead=sscanf(line, "%s %s %s %s %s %s %s %s ",newName,sa[0],sa[1],sa[2],sa[3],sa[4],sa[5],sa[6]); break;
	case 8: numRead=sscanf(line, "%s %s %s %s %s %s %s %s %s ",newName,sa[0],sa[1],sa[2],sa[3],sa[4],sa[5],sa[6],sa[7]); break;
	default: Fatal("hmm, unknown value of k %d", k); break;
	}
	if(numRead != k+1) Fatal("for k=%d we expect %d columns but got %d", k,k+1,numRead);
	if(strcmp(name, newName) != 0)
	{
	    strcpy(name, newName);
	    for(i=0;i<k;i++) strcpy(Sarray[i], sa[i]);
	}
	for(i=0; i<k-1; i++)for(j=i+1;j<k;j++)
	{
	    assert(i!=j);
	    int e1 = GraphAreConnected(G, Gni(G, sa[i]),      Gni(G,     sa[j]));
	    int e2 = GraphAreConnected(G, Gni(G, Sarray[i]),  Gni(G, Sarray[j]));
	    if(e1 == e2) /* do nothing */ ;
	    else {
		printf("line %d: edge(%s,%s)=%d, edge(%s,%s)=%d\t", lines,
		    sa[i], sa[j], e1, Sarray[i], Sarray[j], e2); fflush(stdout);
		fputs(line, stdout); fflush(stdout);
		//assert(false);
	    }
	}
    }
    if(n != lines) Fatal("blant-sanity: expected %d lines but got %d\n",n,lines);
    return 0;
}

