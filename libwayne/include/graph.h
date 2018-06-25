#ifndef _GRAPH_H
#define _GRAPH_H

#include "misc.h"
#include "sets.h"
#include "combin.h"
#include <stdio.h>

#define SORT_NEIGHBORS 0 // Thought this might speed things up but it appears not to.

#define SUPPORT_NODE_NAMES 1 // allows node names to be strings but slows things down alot for big graphs
#if SUPPORT_NODE_NAMES
#include "bintree.h"
#endif

/* Constructs for simple graphs, no self-loops: edge (i,i) never exists.
*/

typedef struct _Graph {
    /* vertices numbered 0..n-1 inclusive */
    int n;
    SET **A;   /* Adjacency Matrix, as a dynamically allocated array[G->n] of SETs */
    Boolean sparse; // true=only neighbors and degree, no matrix; false=only matrix + degree, no neighbors
    int *degree;   /* degree of each v[i] == cardinality of A[i] == length of neighbor array */
    int **neighbor; /* adjacency list: possibly sorted list of neighbors, sorted if SORTED below is true. */
#if SORT_NEIGHBORS
    SET *sorted; // Boolean array: when sparse, is the neighbor list of node[i] sorted or not?
#endif
    int maxEdges, numEdges, *edgeList; /* UNSORTED list of all edges in the graph, edgeList[0,..2*numEdges] */
#if SUPPORT_NODE_NAMES
    BINTREE *nameDict;	// string to int map
    char **name;	// int to string map (inverse of the above)
#endif
} GRAPH;

GRAPH *GraphAlloc(unsigned int n, Boolean sparse);
void GraphFree(GRAPH *G);
GRAPH *GraphEdgesAllDelete(GRAPH *G);
GRAPH *GraphConnect(GRAPH *G, int i, int j);
GRAPH *GraphDisconnect(GRAPH *G, int i, int j);
GRAPH *GraphComplement(GRAPH *Gbar, GRAPH *G);
GRAPH *GraphUnion(GRAPH *destination, GRAPH *G1, GRAPH *G2);
#define GraphDegree(G,v) ((G)->degree[v])
GRAPH *GraphCopy(GRAPH *Gc, GRAPH *G); // If Gc == NULL, create duplicate.  Otherwise just copy G's info into Gc.


/* Returns number of nodes in the the distance-d neighborhood, including seed.
 * nodeArray[0] always= seed, distArray[seed] always= 0.
 * distance = 0 implies no BFS, seed as only member of nodeArray, and return value 1.
 * Use distance >= G->n to ensure BFS will go out as far as possible.
 * Both arrays should have at least G->n elements; caller is reponsible for allocating them.
 * For each nodeArray[i], distArray[nodeArray[i]] is the distance from seed.
 * (Note that distArray[j] has an entry for EVERY node in the graph; only those
 * with indices j \in nodeArray have meaning.)
 */
int GraphBFS(GRAPH *G, int seed, int distance, int *nodeArray, int *distArray);

/* Uses DFS on G starting at node v to see if the current connected component has at least k nodes.
** More efficient than a full DFS. */
Boolean GraphCCatLeastK(GRAPH *G, int v, int k);
Boolean _GraphCCatLeastKHelper(GRAPH *G, SET* visited, int v, int *k);

/* Full DFS on whatever connected component v is in.  On top-level call, you should set (*pn)=0.
** We will populate Varray with the elements and also set their visited state to true.
** The Varray will be nuked starting at position (*pn), and visited does *not* need to be clear,
** so you can call this function multiple times for each connected component and in the end all
** the visited values should be true, and Varray will have all the elements, ordered by which
** connected component they are in.  We return *pn at the end.
*/
int GraphVisitCC(GRAPH *G, unsigned int v, SET *visited, unsigned int *Varray, int *pn);


/*
** GraphInduced_NoVertexDelete doesn't delete any vertices, it only deletes
** edges whose ends don't both appear in V.  GraphInduced builds an entirely
** new graph in which also the vertices not in V are deleted; vertices
** are renumbered but ordering is conserved.
*/
GRAPH *GraphInduced_NoVertexDelete(GRAPH *Gi, GRAPH *G, SET *V);
GRAPH *GraphInduced(GRAPH *Gi, GRAPH *G, SET *V);

void GraphPrintAdjMatrix(FILE *fp, GRAPH *G);
GRAPH *GraphReadAdjMatrix(FILE *fp, Boolean sparse);
void GraphPrintAdjList(FILE *fp, GRAPH *G);
GRAPH *GraphReadAdjList(FILE *fp, Boolean sparse);
GRAPH *GraphFromEdgeList(int numNodes, int numEdges, int *pairs, Boolean sparse);
GRAPH *GraphReadEdgeList(FILE *fp, Boolean sparse);
void GraphPrintConnections(FILE *fp, GRAPH *G);
GRAPH *GraphReadConnections(FILE *fp, Boolean sparse);
int GraphNumEdges(GRAPH *G); // total number of edges, just the sum of the degrees / 2.

// Only enable this macro if you're not allowing sparse graphs
//#define GraphAreConnected(G,i,j) SetIn((G)->A[i],(j))
#ifndef GraphAreConnected
Boolean GraphAreConnected(GRAPH *G, int i, int j);
#endif

/*
** The following subroutines should be used with caution, because they take
** exponential time.  They look for cliques and independent sets of size n,
** respectively.  They return 0 if the requested property doesn't exist, or
** else a SET containing the nodes with the property requested.  Calling
** the appropriate GraphNext* function returns the next in the list of
** Kn's or In's. 0 is returned when none are left.
** The underlying graph G should NOT be changed while these routines are
** "in motion".
*/
typedef struct _clique {
    GRAPH *G;
    SET *set;
    unsigned cliqueSize, *combArray, *inducedArray;
    COMBIN* combin;
} CLIQUE;
CLIQUE *GraphKnFirst(GRAPH *G, int n);
Boolean GraphKnContains(GRAPH *G, int n);
SET *GraphKnNext(CLIQUE*);
CLIQUE *GraphInFirst(GRAPH *G, int n);
Boolean GraphInContains(GRAPH *G, int n);
#define GraphInNext GraphKnNext
void GraphCliqueFree(CLIQUE *);

/*
** Isomorphism algorithm for general graphs G1, G2.
** This algorithm uses some simple tests in an attempt to avoid the
** exponential algorithm.  It checks the number of nodes, the number of
** times each degree appears, and then kranks up the exponential time
** algorithm.
*/
Boolean GraphsIsomorphic(int *perm, GRAPH *G1, GRAPH *G2);

#endif /* _GRAPH_H */
