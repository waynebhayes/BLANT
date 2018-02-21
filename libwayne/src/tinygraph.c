#include "misc.h"
#include "sets.h"
#include "tinygraph.h"
#include "queue.h"
#include <assert.h>

/*************************************************************************
**
**                            The Basics
**
*************************************************************************/

TINY_GRAPH *TinyGraphAlloc(unsigned int n)
{
    static Boolean startup = 1;
    TINY_GRAPH *G = Calloc(1, sizeof(TINY_GRAPH));
    assert(n <= MAX_TSET);
    if(startup)
    {
	startup = 0;
	SetStartup();
    }
    G->n = n;
    return G;
}

TINY_GRAPH *TinyGraphConnect(TINY_GRAPH *G, int i, int j)
{   
    if(TinyGraphAreConnected(G, i, j))
	return G;
    TSetAdd(G->A[i], j);
    TSetAdd(G->A[j], i);
    ++G->degree[i];
    ++G->degree[j];
    return G;
}

TINY_GRAPH *TinyGraphEdgesAllDelete(TINY_GRAPH *G)
{
    int i;
    for(i=0; i < MAX_TSET; i++)
    {
	G->degree[i] = 0;
	TSetEmpty(G->A[i]);
    }
    return G;
}

TINY_GRAPH *TinyGraphDisconnect(TINY_GRAPH *G, int i, int j)
{   
    if(!TinyGraphAreConnected(G, i, j))
	return G;
    TSetDelete(G->A[i], j);
    TSetDelete(G->A[j], i);
    --G->degree[i];
    --G->degree[j];
    return G;
}

#ifndef TinyGraphAreConnected
Boolean TinyGraphAreConnected(TINY_GRAPH *G, int i, int j)
{
    if(TSetIn(G->A[i],j))
    {
	assert(TSetIn(G->A[j],i));
	return true;
    }
    else
	return false;
}
#endif

void TinyGraphPrintAdjMatrix(FILE *fp, TINY_GRAPH *G)
{
    int i, j;
    for(i=0; i<G->n; i++)
    {
	for(j=0; j<G->n - 1; j++)
	    fprintf(fp, "%d ", !!TinyGraphAreConnected(G,i,j));
	fprintf(fp, "%d\n", !!TinyGraphAreConnected(G,i,j));
    }
}


TINY_GRAPH *TinyGraphReadAdjMatrix(FILE *fp)
{
    int i,j,n;
    TINY_GRAPH *G;
    fscanf(fp, "%d", &n);
    G = TinyGraphAlloc(n);
    for(i=0; i<n; i++) for(j=0; j<n; j++)
    {
	int connected;
	if(fscanf(fp, "%d", &connected) != 1)
	    Fatal("too few vertices listed in file");
	if(connected)
	    TinyGraphConnect(G,i,j);
    }
    return G;
}


TINY_GRAPH *TinyGraphComplement(TINY_GRAPH *Gbar, TINY_GRAPH *G)
{
    int i;
    TSET full = 0;
    if(!Gbar)
	Gbar = TinyGraphAlloc(G->n);

    for(i=0; i < G->n; i++)
	TSetAdd(full, i);

    Gbar->n = G->n;

    /* A tad tricky: we don't want to include self-loops, so we add them,
    ** and then XOR with the full row.
    */
    for(i=0; i < G->n; i++)
    {
	Gbar->A[i] = (G->A[i] | (TSET1 << i)) ^ full;
	Gbar->degree[i] = G->n - 1 - G->degree[i];
    }
    return Gbar;
}


TINY_GRAPH *TinyGraphUnion(TINY_GRAPH *dest, TINY_GRAPH *G1, TINY_GRAPH *G2)
{
    int i;

    if(G1->n != G2->n)
	return NULL;
    
    if(dest)
	dest->n = G1->n;
    else
	dest = TinyGraphAlloc(G1->n);

    for(i=0; i < G1->n; i++)
    {
	dest->A[i] = G1->A[i] | G2->A[i];
	dest->degree[i] = TSetCardinality(dest->A[i]);
    }

    return dest;
}


TINY_GRAPH *TinyGraphCopy(TINY_GRAPH *dest, TINY_GRAPH *G1)
{
    int i;

    if(dest)
	dest->n = G1->n;
    else
	dest = TinyGraphAlloc(G1->n);

    for(i=0; i < G1->n; i++)
    {
	dest->A[i] = G1->A[i];
	dest->degree[i] = G1->degree[i];
    }

    return dest;
}


int TinyGraphBFS(TINY_GRAPH *G, int root, int distance, int *nodeArray, int *distArray)
{
    QUEUE *BFSQ;
    int i, count = 0;

    assert(0 <= root && root < G->n);
    assert(distance >= 0);
    assert(nodeArray != NULL);
    assert(distArray != NULL);

    if(distance == 0) /* We could let the rest of the routine run, but why bother? */
    {
	nodeArray[0] = root;
	distArray[root] = 0;
	return 1;
    }

    for(i=0; i<G->n; i++)
	nodeArray[i] = distArray[i] = -1;

    distArray[root] = 0;
    BFSQ = QueueAlloc(G->n);
    QueuePut(BFSQ, (foint)root);
    while(QueueSize(BFSQ) > 0)
    {
	int v = QueueGet(BFSQ).i;

	/* At this point, distArray[v] should be assigned (when v was appended
	 * onto the queue), but v hasn't been "visited" or "counted" yet.
	 */

	assert(0 <= v && v < G->n);
	assert(0 <= distArray[v] && distArray[v] < G->n);

	assert(nodeArray[count] == -1);
	nodeArray[count] = v;
	count++;

	if(distArray[v] < distance) /* v's neighbors will be within BFS distance */
	{
	    int neighbor[MAX_TSET];
	    int j, numNeighbors = TSetToArray(neighbor, G->A[v]); /* This is the slow part, O(n) */
	    for(j=0; j < numNeighbors; j++)
		if(distArray[neighbor[j]] == -1) /* some of the neighbors might have already been visited */
		{
		    distArray[neighbor[j]] = distArray[v] + 1;
		    QueuePut(BFSQ, (foint)neighbor[j]);
		}
	}
    }
    QueueFree(BFSQ);
    return count;
}

TINY_GRAPH *TinyGraphInduced(TINY_GRAPH *Gv, TINY_GRAPH *G, TSET V)
{
    unsigned int array[MAX_TSET], nV = TSetToArray(array, V), i, j;
    TINY_GRAPH GGv;
    if(Gv)
    {
	if(Gv == G) /* be careful, destination is same as source */
	    Gv = &GGv;
	Gv->n = nV;
	TinyGraphEdgesAllDelete(Gv);
    }
    else
	Gv = TinyGraphAlloc(nV);

    for(i=0; i < nV; i++) for(j=i+1; j < nV; j++)
	if(TinyGraphAreConnected(G, array[i], array[j]))
	    TinyGraphConnect(Gv, i, j);

    if(Gv == &GGv)
	*(Gv = G) = GGv;
    return Gv;
}


TINY_GRAPH *TinyGraphInduced_NoVertexDelete(TINY_GRAPH *Gv, TINY_GRAPH *G, TSET V)
{
    unsigned int array[MAX_TSET], nV = TSetToArray(array, V), i, j;
    TINY_GRAPH GGv;
    if(Gv)
    {
	if(Gv == G) /* be careful, destination is same as source */
	    Gv = &GGv;
	Gv->n = G->n;
	TinyGraphEdgesAllDelete(Gv);
    }
    else
	Gv = TinyGraphAlloc(G->n);

    for(i=0; i < nV; i++) for(j=i+1; j < nV; j++)
	if(TinyGraphAreConnected(G, array[i], array[j]))
	    TinyGraphConnect(Gv, array[i], array[j]);
    if(Gv == &GGv)
	*(Gv = G) = GGv;
    return Gv;
}

/*
** A reasonably fast search for triangles.  If there exists a triangle,
** it'll return the TSET representation of it, otherwise return 0.
*/
#if 0
static int TinyGraphContainsK3(TINY_GRAPH *G)
{
    int i,j;
    for(i=0; i < G->n; i++)
	for(j=i+1; j < G->n; j++)
	    if(TinyGraphAreConnected(G,i,j) && TinyGraphConnectingCausesK3(G,i,j))
	    {
		int k;
		for(k=0; k < G->n; k++)
		    if(TSetIn(TSetIntersect(G->A[i], G->A[j]), k))
			break;
		return (TSET1<<i)|(TSET1<<j)|(TSET1<<k);
	    }
    return 0;
}
#endif


/**************************************************************************
**
**  Graph Isomorphism
**
**************************************************************************/

static TINY_GRAPH *isoG1, *isoG2;

static int _permutationIdentical(int n, int perm[n])
{
    int i, j;
    for(i=0; i<n; i++)
	if(isoG1->degree[i] != isoG2->degree[perm[i]])
	    return 0;

    for(i=0; i<n; i++) for(j=i+1; j<n; j++)
	/* The !GraphAreConnected is just to turn a bitstring into a boolean */
	if(!TinyGraphAreConnected(isoG1, i,j) !=
	    !TinyGraphAreConnected(isoG2, perm[i], perm[j]))
	    return 0;   /* non-isomorphic */
    return 1;   /* isomorphic! */
}

/* If returns true, then populates the perm array with the permutation proving isomorphism
** You need to provide the perm array long enough (of size G1->n).
*/
Boolean TinyGraphsIsomorphic(int *perm, TINY_GRAPH *G1, TINY_GRAPH *G2)
{
    int i, n = G1->n, degreeCount1[n], degreeCount2[n];
    TSET degreeOnce;

    /*
    ** First some simple tests.
    */
    if(G1->n != G2->n)
	return false;
    
    if(n < 2)
	return true;

    /*
    ** Ensure each degree occurs the same number of times in each.
    */
    for(i=0; i<n; i++)
	degreeCount1[i] = degreeCount2[i] = 0;
    for(i=0; i<n; i++)
    {
	++degreeCount1[G1->degree[i]];
	++degreeCount2[G2->degree[i]];
    }
    for(i=0; i<n; i++)
	if(degreeCount1[i] != degreeCount2[i])
	    return false;

    /*
    ** Let degree d appear only once.  Then there is exactly one vertex
    ** v1 in G1 with degree d, and exactly one vertex v2 in G2 with degree d.
    ** G1 and G2 are isomorphic only if the neighborhood of v1 is isomorphic
    ** to the neighborhood of v2, and also if (G1-v1) is isomorphic to (G2-v2).
    */
    TSetEmpty(degreeOnce);
    for(i=0; i<n; i++)
	if(degreeCount1[i] == 1)
	    TSetAdd(degreeOnce, i);

    for(i=0; i<n; i++)
    {
	/* Find out if the degree of vertex i in G1 appears only once */
	if(TSetIn(degreeOnce, G1->degree[i]))
	{
	    int j, degree = G1->degree[i];
	    TINY_GRAPH neighG1i, neighG2j, restG1i, restG2j;

	    /* find the (unique) vertex in G2 that has the same degree */
	    for(j=0; j < n; j++)
		if(G2->degree[j] == degree)
		    break;
	    assert(j < n);

	    TinyGraphInduced(&neighG1i, G1, G1->A[i]);
	    TinyGraphInduced(&neighG2j, G2, G2->A[j]);

	    /*
	    ** Note: this recursion works only as long as
	    ** _permutationIdentical doesn't call GraphsIsomorphic.
	    ** (if it does, isoG1 and isoG2 get messed up).
	    */

	    if(!TinyGraphsIsomorphic(perm, &neighG1i, &neighG2j))
		return false;

	    /* Now ask if the remainder of the graphs are isomorphic */
	    TSET all = ((unsigned char)-1) >> (8-G1->n); 
	    TinyGraphInduced(&restG1i, G1, TSetDelete(all, i));
	    all = ((unsigned char)-1) >> (8-G1->n); 
	    TinyGraphInduced(&restG2j, G2, TSetDelete(all, j));
	    if(!TinyGraphsIsomorphic(perm, &restG1i, &restG2j))
		return false;
	}
    }

    /*
    ** Oh well, fire up the exponential search.
    ** CombinAllPermutations will return 0 iff all permutations were
    ** tried; the function _permutationIdentical should return non-zero
    ** when it finds an identical permutation, and that non-zero value
    ** will be returned here, indicating an identical permutation was
    ** found, ie, that the graphs are isomorphic.
    */
    isoG1 = G1; isoG2 = G2;
    return !!CombinAllPermutations(n, perm, _permutationIdentical);
}
