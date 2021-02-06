#ifndef _TINYGRAPH_H
#define _TINYGRAPH_H

#include "misc.h"
#include "sets.h"
#include "combin.h"
#include <stdio.h>

/* Constructs for simple graphs, no self-loops: edge (i,i) never exists.
*/

typedef struct _tinyGraph {
    /* vertices numbered 0..n-1 inclusive; n must be <= MAX_TSET */
    char n, degree[MAX_TSET];   /* degree of each v[i] == cardinality of A[i] */
    TSET A[MAX_TSET];   /* Adjacency Matrix */
} TINY_GRAPH;

TINY_GRAPH *TinyGraphAlloc(unsigned int n);
#define TinyGraphFree free
TINY_GRAPH *TinyGraphEdgesAllDelete(TINY_GRAPH *G);
TINY_GRAPH *TinyGraphCopy(TINY_GRAPH *G, TINY_GRAPH *H); // G = H
TINY_GRAPH *TinyGraphConnect(TINY_GRAPH *G, int i, int j);
TINY_GRAPH *TinyGraphDisconnect(TINY_GRAPH *G, int i, int j);
TINY_GRAPH *TinyGraphComplement(TINY_GRAPH *Gbar, TINY_GRAPH *G);
TINY_GRAPH *TinyGraphUnion(TINY_GRAPH *destination, TINY_GRAPH *G1, TINY_GRAPH *G2);
int TinyGraphNumEdges(TINY_GRAPH *G); // total number of edges, just the sum of the degrees / 2.
#define TinyGraphDegree(G,v) ((G)->degree[v])

/* Returns number of nodes in the the distance-d neighborhood, including seed.
 * nodeArray[0] always= seed, distArray[seed] always= 0.
 * distance = 0 implies no BFS, seed as only member of nodeArray, and return value 1.
 * Use distance >= G->n to ensure BFS will go out as far as possible.
 * Both arrays should have at least G->n elements; caller is reponsible for allocating them.
 * For each nodeArray[i], distArray[nodeArray[i]] is the distance from seed.
 * (Note that distArray[j] has an entry for EVERY node in the graph; only those
 * with indices j \in nodeArray have meaning.)
 * Also, worst-case runtime is O(n^2)... very bad, yes.. :-(
 */
int TinyGraphBFS(TINY_GRAPH *G, int seed, int distance, int *nodeArray, int *distArray);
Boolean TinyGraphDFSConnected(TINY_GRAPH *G, int seed);
void TinyGraphDFSConnectedHelper(TINY_GRAPH *G, int seed, TSET* visited);
/*
** TinyGraphInduced_NoVertexDelete doesn't delete any vertices, it only deletes
** edges whose ends don't both appear in V.  TinyGraphInduced builds an entirely
** new graph in which also the vertices not in V are deleted; vertices
** are renumbered but ordering is conserved.
*/
TINY_GRAPH *TinyGraphInduced_NoVertexDelete(TINY_GRAPH *Gi, TINY_GRAPH *G, TSET V);
TINY_GRAPH *TinyGraphInduced(TINY_GRAPH *Gi, TINY_GRAPH *G, TSET V);

void TinyGraphPrintAdjMatrix(FILE *fp, TINY_GRAPH *G);
TINY_GRAPH *TinyGraphReadAdjMatrix(FILE *fp);

#define TinyGraphAreConnected(G,i,j) TSetIn((G)->A[i],(j))
#ifndef TinyGraphAreConnected
Boolean TinyGraphAreConnected(TINY_GRAPH *G, int i, int j);
#endif

/*
** Isomorphism algorithm for general graphs G1, G2.
** This algorithm uses some simple tests in an attempt to avoid the
** exponential algorithm.  It checks the number of nodes, the number of
** times each degree appears, and then kranks up the exponential time
** algorithm.
*/
Boolean TinyGraphsIsomorphic(int *perm, TINY_GRAPH *G1, TINY_GRAPH *G2);

#endif /* _TINYGRAPH_H */
