#include "misc.h"
#include "sets.h"
#include "graph.h"
#include "queue.h"
#include "rand48.h"
#include "Oalloc.h"
#include <assert.h>
#include <string.h>

#define MIN_EDGELIST 1024

/*************************************************************************
**
**                            The Basics
**
*************************************************************************/

GRAPH *GraphAlloc(unsigned int n, Boolean sparse, Boolean supportNodeNames)
{
    static Boolean needStartup = 1;
    int i;
    GRAPH *G = Calloc(1, sizeof(GRAPH));
    if(needStartup)
    {
	needStartup = 0;
	SetStartup();
    }
    G->sparse = sparse;
    G->n = n;
    G->A = NULL;
    G->degree = Calloc(n, sizeof(G->degree[0]));
    G->maxEdges = MIN_EDGELIST;
    if(sparse>=true) // >=true means "both"
    {
	G->edgeList = Malloc(2*G->maxEdges*sizeof(int));
	G->numEdges = 0;
	G->neighbor = Calloc(n, sizeof(G->neighbor[0]));
#if SORT_NEIGHBORS
	G->sorted = SetAlloc(G->n);
#endif
    }
    if(sparse==both || sparse == false)
    {
	G->A = Calloc(n, sizeof(G->A[0]));
	for(i=0; i<n; i++)
	    G->A[i] = SetAlloc(n);
    }
    G->supportNodeNames = supportNodeNames;
    G->nameDict = NULL;
    return G;
}

void GraphFree(GRAPH *G)
{
    int i;
    for(i=0; i<G->n; i++)
    {
	if(G->sparse>=true) Free(G->neighbor[i]);
	if(!G->sparse||G->sparse==both)
	    SetFree(G->A[i]);
    }
    Free(G->degree);
    if(G->sparse >= true)
    {
	Free(G->edgeList);
	Free(G->neighbor);
    }
    if(!G->sparse||G->sparse==both)
	Free(G->A);
    if(G->name) {
	for(i=0;i<G->n;i++) Free(G->name[i]);
	Free(G->name);
    }
    if(G->nameDict) BinTreeFree(G->nameDict);
    Free(G);
}

/* If Gc == NULL, create duplicate.  Otherwise just copy G's info into Gc. */
GRAPH *GraphCopy(GRAPH *Gc, GRAPH *G)
{
    int i;
    if(Gc == G) return G;
    if(!Gc)
	Gc = GraphAlloc(G->n, G->sparse, G->supportNodeNames);

    Gc->sparse = G->sparse;
    if(G->n > Gc->n)
    {
	if(G->sparse >= true)
	{
	    Gc->A = NULL;
	    Gc->neighbor = Realloc(Gc->neighbor, G->n * sizeof(Gc->neighbor[0]));
	    Gc->numEdges = G->numEdges;
	    Gc->maxEdges = G->maxEdges;
	    Gc->edgeList = Realloc(Gc->edgeList, 2*G->maxEdges*sizeof(int));
	    for(i=0;i < G->numEdges; i++)
	    {
		Gc->edgeList[2*i] = G->edgeList[2*i];
		Gc->edgeList[2*i+1] = G->edgeList[2*i+1];
	    }
	}
	if(!G->sparse || G->sparse==both)
	    Gc->A = Realloc(Gc->A, G->n * sizeof(Gc->A[0]));
	Gc->degree = Realloc(Gc->degree, G->n * sizeof(Gc->degree[0]));
	/* reallaoc doesn't zero out new entries, damn it */
	for(i=Gc->n; i < G->n; i++)
	{
	    if(Gc->sparse >= true)
		Gc->neighbor[i] = NULL;
	    if(!Gc->sparse || Gc->sparse==both)
		Gc->A[i] = NULL;
	}
    }

    for(i=0; i < G->n; i++)
    {
	if((!Gc->sparse || Gc->sparse==both) && Gc->A[i] && Gc->n != G->n)
	{
	    SetFree(Gc->A[i]);
	    Gc->A[i] = NULL;
	}
	if(!Gc->sparse||Gc->sparse==both) Gc->A[i] = SetCopy(Gc->A[i], G->A[i]);
	else if(G->degree[i] > Gc->degree[i])
	    Gc->neighbor[i] = Realloc(Gc->neighbor[i], G->degree[i] * sizeof(Gc->neighbor[i][0]));
	if(Gc->sparse>=true) memmove(Gc->neighbor[i], G->neighbor[i], G->degree[i] * sizeof(G->neighbor[i][0]));
	Gc->degree[i] = G->degree[i];
    }

    Gc->n = G->n;
    return Gc;
}

#if SORT_NEIGHBORS
// Used when qsort'ing the neighbors when graph is sparse.
static int IntCmp(const void *a, const void *b)
{
    const int *i = (const int*)a, *j = (const int*)b;
    return (*i)-(*j);
}

static GRAPH *GraphSort(GRAPH *G)
{
    if(G->sparse>=true)
    {
	int v;
	for(v=0; v<G->n; v++) if(!SetIn(G->sorted, v))
	{
	    qsort(G->neighbor[v], G->degree[v], sizeof(G->degree[0]), IntCmp);
	    SetAdd(G->sorted, v);
	}
    }
    return G;
}
#else
#define GraphSort(x)
#endif

GRAPH *GraphConnect(GRAPH *G, int i, int j)
{
    assert(0 <= i && i < G->n && 0 <= j && j < G->n && i!=j);
    if(GraphAreConnected(G, i, j))
    {
	assert(GraphAreConnected(G, j, i));
	return G;
    }
    if(G->sparse>=true)
    {
	G->neighbor[i] = Realloc(G->neighbor[i], (G->degree[i]+1)*sizeof(int));
	G->neighbor[j] = Realloc(G->neighbor[j], (G->degree[j]+1)*sizeof(int));
	assert(G->neighbor[i]);
	assert(G->neighbor[j]);
	G->neighbor[i][G->degree[i]] = j;
	G->neighbor[j][G->degree[j]] = i;
#if SORT_NEIGHBORS
	SetDelete(G->sorted, i);
	SetDelete(G->sorted, j);
#endif
	assert(G->numEdges <= G->maxEdges);
	if(G->numEdges == G->maxEdges)
	{
	    G->maxEdges *= 2;
	    G->edgeList = Realloc(G->edgeList, 2*G->maxEdges*sizeof(int));
	    assert(G->edgeList);
	}
	G->edgeList[2*G->numEdges] = MIN(i,j);
	G->edgeList[2*G->numEdges+1] = MAX(i,j);
	G->numEdges++;
    }
    if(!G->sparse || G->sparse==both)
    {
	SetAdd(G->A[i], j);
	SetAdd(G->A[j], i);
    }
    ++G->degree[i];
    ++G->degree[j];
    return G;
}

GRAPH *GraphEdgesAllDelete(GRAPH *G)
{
    int i;
    for(i=0; i < G->n; i++)
    {
	if(!G->sparse||G->sparse==both) SetEmpty(G->A[i]);
	G->degree[i] = 0;
#if SORT_NEIGHBORS
	SetDelete(G->sorted, i);
#endif
	/* Don't need to realloc/free neighbors, it'll happen automatically once we start re-adding edges */
    }
    G->numEdges = 0;
    G->maxEdges = MIN_EDGELIST;
    if(G->sparse>=true) G->edgeList = Realloc(G->edgeList, 2*G->maxEdges*sizeof(int));
    return G;
}

GRAPH *GraphDisconnect(GRAPH *G, int i, int j)
{
    int k;
    assert(0 <= i && i < G->n && 0 <= j && j < G->n);
    if(!GraphAreConnected(G, i, j))
	return G;
    --G->degree[i];
    --G->degree[j];

    if(!G->sparse||G->sparse==both)
    {
	SetDelete(G->A[i], j);
	SetDelete(G->A[j], i);
	return G;
    }

    assert(G->sparse>=true);
    if(j < i)
    {
	int tmp = i; i=j; j=tmp;
    }

    Boolean found=false;
    for(k=0; k < G->numEdges; k++)
    {
	if(G->edgeList[2*k] == i && G->edgeList[2*k+1]==j)
	{
	    G->numEdges--;
	    G->edgeList[2*k] = G->edgeList[2*G->numEdges];
	    G->edgeList[2*k+1] = G->edgeList[2*G->numEdges+1];
	    found=true;
	    break;
	}
    }
    assert(found);

    /* now find and delete each other's neighbors */
    k=0;
    while(G->neighbor[i][k] != j)
	k++;
    assert(k <= G->degree[i] && G->neighbor[i][k] == j); /* this is the new degree, so using "<=" is correct */
    G->neighbor[i][k] = G->neighbor[i][G->degree[i]];

    k=0;
    while(G->neighbor[j][k] != i)
	k++;
    assert(k <= G->degree[j] && G->neighbor[j][k] == i);
    G->neighbor[j][k] = G->neighbor[j][G->degree[j]];
#if SORT_NEIGHBORS
    SetDelete(G->sorted, i);
    SetDelete(G->sorted, j);
#endif
    return G;
}

#ifndef GraphAreConnected
Boolean GraphAreConnected(GRAPH *G, int i, int j)
{
    assert(0 <= i && i < G->n && 0 <= j && j < G->n);
    if(G->sparse>=true)
    {
#if SORT_NEIGHBORS
	if(SetIn(G->sorted, i))
	    return !!bsearch(&j, G->neighbor[i], G->degree[i], sizeof(G->neighbor[0]), IntCmp);
	else if(SetIn(G->sorted, j))
	    return !!bsearch(&i, G->neighbor[j], G->degree[j], sizeof(G->neighbor[0]), IntCmp);
	else
#endif
	{
	    int k, n, *neighbors, me, other;
	    // Check through the shorter list
	    if(G->degree[i] < G->degree[j])
	    {
		me = i; other = j;
	    }
	    else
	    {
		me = j; other = i;
	    }
	    n = G->degree[me];
	    neighbors = G->neighbor[me];
	    for(k=0; k<n; k++)
		if(neighbors[k] == other)
		    return true;
	    return false;
	}
    }
    if(!G->sparse||G->sparse==both)
    {
	if(SetIn(G->A[i],j))
	{
	    if(!SetIn(G->A[j],i))
		Fatal("SetIn(%d,%d)=%ld, SetIn(%d,%d)=%ld\n", i,j,SetIn(G->A[i],j), j,i,SetIn(G->A[j],i));
	    return true;
	}
	else
	    return false;
    }
    return false;
}
#endif


int GraphNumCommonNeighbors(GRAPH *G, int i, int j)
{
    assert(0 <= i && i < G->n && 0 <= j && j < G->n);
    int numCommon1 = 0, numCommon2 = 0;
    if(G->sparse>=true) // loop version is MUCH faster, so use it if possible
    {
	int k, n, *neighbors, me, other;
	// Check through the shorter list
	if(G->degree[i] < G->degree[j])
	{
	    me = i; other = j;
	}
	else
	{
	    me = j; other = i;
	}
	n = G->degree[me];
	neighbors = G->neighbor[me];
	for(k=0; k<n; k++)
	    if(neighbors[k] != other && GraphAreConnected(G, neighbors[k], other))
		++numCommon1;
    }
    else
    {
	assert(!G->sparse||G->sparse==both);
	static SET *intersect;
	if(!intersect) intersect = SetAlloc(G->n);
	SetReset(intersect);
	SetIntersect(intersect, G->A[i], G->A[j]);
#if PARANOID_ASSERTS
	assert(!SetIn(intersect,i) && !SetIn(intersect,j));
#endif
	numCommon2 = SetCardinality(intersect); // no self-loops, so this works even if edge (i,j) exists (assertion above)
	if(numCommon1) {
	    assert(numCommon1 == numCommon2);
	    numCommon1=0;
	}
    }
    assert(numCommon1 == 0 || numCommon2 == 0);
    return numCommon1 + numCommon2;
}

int GraphNumEdges(GRAPH *G)
{
    int total=0, i;
    for(i=0; i<G->n; i++)
	total += G->degree[i];
    assert(total % 2 == 0); // should be divisible by 2
    assert(G->numEdges == total/2);
    return G->numEdges;
}

void GraphPrintAdjMatrix(FILE *fp, GRAPH *G)
{
    int i, j;
    fprintf(fp, "%d\n", G->n);
    for(i=0; i<G->n; i++)
    {
	for(j=0; j<G->n - 1; j++)
	    fprintf(fp, "%d ", !!GraphAreConnected(G,i,j));
	fprintf(fp, "%d\n", !!GraphAreConnected(G,i,j));
    }
}


GRAPH *GraphReadAdjMatrix(FILE *fp, Boolean sparse)
{
    int i,j,n;
    GRAPH *G;
    if(fscanf(fp, "%d", &n) != 1)
	Fatal("GraphReadAdjMatrix: reading 'n' failed");
    assert(n >= 0);
    G = GraphAlloc(n, sparse, false); // no SUPPORT_NODE_NAMES at the moment
    for(i=0; i<n; i++) for(j=0; j<n; j++)
    {
	int connected;
	if(fscanf(fp, "%d", &connected) != 1)
	    Fatal("GraphReadAdjMatrix: reading entry(%d,%d) failed", i, j);
	if(connected)
	    GraphConnect(G,i,j);
    }
    GraphSort(G);
    return G;
}


void GraphPrintAdjList(FILE *fp, GRAPH *G)
{
    int i, j;
    assert(G->sparse>=true);
    fprintf(fp, "%d\n", G->n);
    for(i=0; i<G->n; i++)
    {
	fprintf(fp, "%d ", G->degree[i]);
	for(j=0; j<G->degree[i] - 1; j++)
	    fprintf(fp, "%d ", G->neighbor[i][j]);
	fprintf(fp, "%d\n", G->neighbor[i][j]);
    }
}


GRAPH *GraphReadAdjList(FILE *fp, Boolean sparse)
{
    GRAPH *G;
    int n, i, j, d;
    if(fscanf(fp, "%d", &n) != 1)
	Fatal("GraphReadAdjList: failed to read 'n'");
    assert(n >= 0);
    G = GraphAlloc(n, sparse, false); // no SUPPORT_NODE_NAMES at the moment
    for(i=0; i<n; i++)
    {
	if(fscanf(fp, "%d", &d) != 1)
	    Fatal("node %d: expecting degree, but couldn't find an integer", i);
	for(j=0; j<d; j++)
	{
	    int neigh;
	    if(fscanf(fp, "%d", &neigh) != 1)
		Fatal("node %d, degree %d, ran out of integers on neighbor #%d", i, d, j);
	    else
	    {
		assert(0 <= neigh && neigh < n);
		assert(neigh != i);
		GraphConnect(G, i, neigh);
	    }
	}
    }
    GraphSort(G);
    return G;
}

GRAPH *GraphFromEdgeList(int n, int m, int *pairs, Boolean sparse)
{
    int i;
    GRAPH *G = GraphAlloc(n, sparse, false); // will set names later
    assert(n == G->n);
    assert(G->degree);
    if(sparse>=true)
	for(i=0;i<n;i++)
	    assert(!G->neighbor[i]);
    for(i=0; i<m; i++)
	GraphConnect(G, pairs[2*i], pairs[2*i+1]);
    if(sparse>=true)
    {
	assert(G->neighbor);
	GraphSort(G);
    }
    return G;
}

// Returns a *constant* string; you need to dup it if you want to keep it
char *HashString(char *s)
{
    static char *hash;
    static int buflen;
    int n = strlen(s);
    if(buflen < n+1){
	buflen = 2*n+1;
	hash = realloc(hash, buflen); // realloc accepts NULL
	assert(hash);
    }
    return hash;
}

GRAPH *GraphReadEdgeList(FILE *fp, Boolean sparse, Boolean supportNodeNames)
{
    int numNodes=0;
    int numEdges=0, maxEdges=MIN_EDGELIST; // these will be increased as necessary during reading
    int *pairs = Malloc(2*maxEdges*sizeof(int));

    // SUPPORT_NODE_NAMES
    int maxNodes=MIN_EDGELIST;
    char **names = NULL;
    BINTREE *nameDict = NULL;
    if(supportNodeNames)
    {
	names = Malloc(maxNodes*sizeof(char*));
	nameDict = BinTreeAlloc((pCmpFcn)strcmp, (pFointCopyFcn)strdup, (pFointFreeFcn)free, NULL, NULL);
    }

    char line[BUFSIZ];
    while(fgets(line, sizeof(line), fp))
    {
	// nuke all whitespace, including DOS carriage returns, from the end of the line
	int len = strlen(line);
	while(isspace(line[len-1])) line[--len]='\0';
	int v1, v2;
	assert(numEdges <= maxEdges);
	if(numEdges >= maxEdges)
	{
	    maxEdges *=2;
	    pairs = Realloc(pairs, 2*maxEdges*sizeof(int));
	}
	if(supportNodeNames)
	{
	    assert(numNodes <= maxNodes);
	    if(numNodes+2 >= maxNodes) // -2 for a bit of extra space
	    {
		maxNodes *=2;
		names = Realloc(names, maxNodes*sizeof(char*));
	    }
	    char name1[BUFSIZ], name2[BUFSIZ];
	    if(sscanf(line, "%s%s ", name1, name2) != 2)
		Fatal("GraphReadEdgeList: line %d does not contain 2 strings\n", numEdges);
	    if(strcmp(name1,name2)==0)
		Fatal("GraphReadEdgeList: line %d has self-loop (%s to itself)\n", numEdges,name1);
	    foint f1, f2;
	    if(!BinTreeLookup(nameDict, (foint)name1, &f1))
	    {
		names[numNodes] = strdup(name1);
		f1.i = numNodes++;
		BinTreeInsert(nameDict, (foint)name1, f1);
	    }
	    if(!BinTreeLookup(nameDict, (foint)name2, &f2))
	    {
		names[numNodes] = strdup(name2);
		f2.i = numNodes++;
		BinTreeInsert(nameDict, (foint)name2, f2);
	    }
	    v1 = f1.i; v2 = f2.i;
	}
	else
	{
	    if(sscanf(line, "%d%d ", &v1, &v2) != 2)
		Fatal("GraphReadEdgeList: tried to read pairs number %d but couldn't find 2 ints\n", numEdges);
	    if(v1==v2)
		Fatal("GraphReadEdgeList: line %d has self-loop (%d to itself)\n", numEdges,v1);
	    numNodes = MAX(numNodes, v1);
	    numNodes = MAX(numNodes, v2);
	}
	if(v1 == v2)
	    Fatal("GraphReadEdgeList: line %d: graph cannot cannot have self-loops\n", numEdges);
	pairs[2*numEdges] = v1;
	pairs[2*numEdges+1] = v2;
	if(pairs[2*numEdges] > pairs[2*numEdges+1])
	{
	    int tmp = pairs[2*numEdges];
	    pairs[2*numEdges] = pairs[2*numEdges+1];
	    pairs[2*numEdges+1] = tmp;
	}
	assert(pairs[2*numEdges] < pairs[2*numEdges+1]);
	numEdges++;
    }
    if(supportNodeNames)
    {
	//printf("BINTREE Dictionary Dump\n");
	int i;
	for(i=0; i<numNodes;i++)
	{
	    foint info;
	    assert(BinTreeLookup(nameDict, (foint)names[i], &info));
	    assert(i == info.i);
	    //printf("%d is %s which in turn is %d\n", i, names[i], info.i);
	}
    }
    else
	numNodes++;	// increase it by one since so far it's simply been the biggest integer seen on the input.

    GRAPH *G = GraphFromEdgeList(numNodes, numEdges, pairs, sparse);
    if((G->supportNodeNames = supportNodeNames))
    {
	G->nameDict = nameDict;
	G->name = names;
    }
    Free(pairs);
    assert(G->maxEdges <= maxEdges);
    assert(G->numEdges <= numEdges);
    if(sparse>=true)
    {
	assert(G->neighbor);
	GraphSort(G);
    }
    return G;
}

int GraphNodeName2Int(GRAPH *G, char *name)
{
    foint info;
    assert(BinTreeLookup(G->nameDict, (foint)name, &info));
    return info.i;
}


void GraphPrintConnections(FILE *fp, GRAPH *G)
{
    int i, j;
    assert(G->sparse>=true);
    fprintf(fp, "%d\n", G->n);
    for(i=0; i<G->n; i++) for(j=0; j<G->degree[i]; j++)
	fprintf(fp, "%d %d\n", i, G->neighbor[i][j]);
}


GRAPH *GraphReadConnections(FILE *fp, Boolean sparse)
{
    GRAPH *G;
    int n, i, j, d;
    if(fscanf(fp, "%d", &n) != 1)
	Fatal("GraphReadConnections: failed to read 'n'");
    assert(n >= 0);
    G = GraphAlloc(n, sparse, false);

    while((d=fscanf(fp, "%d %d", &i, &j)) == 2)
    {
	if(i==-1 && j==-1)
	{
	    d=0;
	    break;
	}
	assert(0 <= i && i < G->n);
	assert(0 <= j && j < G->n);
	GraphConnect(G, i, j);
    }
    if(d > 0)
	Fatal("expecting no more integers, but got %d integers", d);
    GraphSort(G);
    return G;
}


GRAPH *GraphComplement(GRAPH *Gbar, GRAPH *G)
{
    int i, j;
    assert(Gbar != G);
    if(!Gbar)
	Gbar = GraphAlloc(G->n, G->sparse, G->supportNodeNames);

    assert(Gbar->n == G->n);

    for(i=0; i < G->n; i++) for(j=i+1; j < G->n; j++)
	if(!GraphAreConnected(G, i, j))
	    GraphConnect(Gbar, i, j);
    GraphSort(Gbar);
    return Gbar;
}


GRAPH *GraphUnion(GRAPH *dest, GRAPH *G1, GRAPH *G2)
{
    int i, j, n = G1->n;

    if(G1->n != G2->n)
	return NULL;

    assert(G1->sparse == G2->sparse);

    if(dest)
	dest->n = n;
    else
	dest = GraphAlloc(n,G1->sparse, G1->supportNodeNames);

    for(i=0; i < n; i++) for(j=i+1; j < n; j++)
	if(GraphAreConnected(G1, i, j) || GraphAreConnected(G2, i, j))
	    GraphConnect(dest, i ,j);
    GraphSort(dest);
    return dest;
}


int GraphBFS(GRAPH *G, int root, int distance, int *nodeArray, int *distArray)
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
	    int j;
	    if(G->sparse>=true)
	    {
		for(j=0; j < G->degree[v]; j++)
		    if(distArray[G->neighbor[v][j]] == -1) /* some of the neighbors might have already been visited */
		    {
			distArray[G->neighbor[v][j]] = distArray[v] + 1;
			QueuePut(BFSQ, (foint)G->neighbor[v][j]);
		    }
	    }
	    if(G->sparse==both || !G->sparse)
	    {
		for(j=0; j < G->n; j++)
		    if(v!=j && GraphAreConnected(G, v, j) && distArray[j] == -1)
		    {
			distArray[j] = distArray[v] + 1;
			QueuePut(BFSQ, (foint)j);
		    }
	    }
	}
    }
    QueueFree(BFSQ);
    return count;
}

Boolean GraphCCatLeastK(GRAPH *G, int v, int k) {
    SET* visited = SetAlloc(G->n);
    Boolean result = _GraphCCatLeastKHelper(G, visited, v, &k);
    SetFree(visited);
    return result;
}

/* visited holds previously visited nodes, v holds the current vertex, k holds the remaining count
** Visit the current node
** if the remaining count reached 0
**      return true
** For each adjacent node
**      if it hasn't been visited
**          recursive call to dfs the node
** return false if the CC wasn't at least k
*/
Boolean _GraphCCatLeastKHelper(GRAPH *G, SET* visited, int v, int *k) {
    SetAdd(visited, v);
    *k -= 1;
    if (*k <= 0) return true;
    int i;
    for (i = 0; i < G->degree[v]; i++) {
        if (!SetIn(visited, G->neighbor[v][i])) {
            Boolean result = _GraphCCatLeastKHelper(G, visited, G->neighbor[v][i], k);
            if (result)
                return result;
        }
    }
    return false;
}

/* At top-level call, set (*pn)=0. The visited array does *not* need to be clear.
** We return the number of elements in Varray.
*/
int GraphVisitCC(GRAPH *G, unsigned int v, SET *visited, unsigned int *Varray, int *pn)
{
    assert(v < visited->n);
    if(!SetIn(visited,v))
    {
	SetAdd(visited, v);
	Varray[(*pn)++] = v;
    	int i;
	for(i=0; i < G->degree[v]; i++)
	    GraphVisitCC(G, G->neighbor[v][i], visited, Varray, pn);
    }
    return *pn;
}

/* doesn't allow Gv == G */
GRAPH *GraphInduced(GRAPH *Gv, GRAPH *G, SET *V)
{
    unsigned array[G->n], nV = SetToArray(array, V), i, j;
    assert(Gv != G);
    if(Gv)
    {
	assert(Gv->n == nV);
	GraphEdgesAllDelete(Gv);
    }
    else
	Gv = GraphAlloc(nV, G->sparse, G->supportNodeNames);
    for(i=0; i < nV; i++) for(j=i+1; j < nV; j++)
	if(GraphAreConnected(G, array[i], array[j]))
	    GraphConnect(Gv, i, j);
    GraphSort(Gv);
    return Gv;
}

/* allows Gv == G */
GRAPH *GraphInduced_NoVertexDelete(GRAPH *Gv, GRAPH *G, SET *V)
{
    unsigned array[G->n], nV = SetToArray(array, V), i, j;
    GRAPH *GGv = GraphAlloc(G->n, G->sparse, G->supportNodeNames);

    for(i=0; i < nV; i++) for(j=i+1; j < nV; j++)
	if(GraphAreConnected(G, array[i], array[j]))
	    GraphConnect(GGv, array[i], array[j]);
    GraphCopy(Gv, GGv);
    GraphFree(GGv);
    GraphSort(Gv);
    return Gv;
}

/*
** A reasonably fast search for triangles.
*/

/*
** This is a helper function for Contains Kn, but you can use it.  It
** tells you if adding this edge will cause a triangle.  Note that
** this function ONLY works if self-loops do not exist!
*/
Boolean GraphConnectingCausesK3(GRAPH *G, int i, int j)
{
    int numIntersect;
    SET *C = SetAlloc(G->n);
    assert(!G->sparse||G->sparse==both);
    SetIntersect(C, G->A[i], G->A[j]);
    numIntersect = SetCardinality(C);
    SetFree(C);
    return numIntersect != 0;
}

Boolean GraphContainsK3(GRAPH *G)
{
    int i,j;
    for(i=0; i < G->n; i++)
	for(j=i+1; j < G->n; j++)
	    if(GraphAreConnected(G,i,j) && GraphConnectingCausesK3(G,i,j))
		return true;
    return false;
}



/**************************************************************************
**
** These are the Clique and Indep set (exponential time) algorithms.
** They work for general graphs.
**
**************************************************************************/

/*
** We check to see if these particular k vertices are a clique by
** bitwise-anding together their adjacency lists (with self-loop added
** to each).  Return false of it's not a clique, otherwise return the set
** of nodes.
** Note that if the graph has changed since the previous call, we may miss some cliques.
*/
SET *Graph_IsCombClique(CLIQUE *c)
{
    SET *intersect;
    int i;

    assert(!c->G->sparse||c->G->sparse==both);
    SetEmpty(c->set);
    for(i=0; i < c->cliqueSize; i++)
	SetAdd(c->set, c->inducedArray[c->combArray[i]]);

    intersect = SetCopy(NULL, c->set);

    for(i=0; i < c->cliqueSize; i++)
    {
	int node = c->inducedArray[c->combArray[i]];
	/* Temporarily add in self-loop for intersection purposes */
	SetAdd(c->G->A[node], node);
	SetIntersect(intersect, intersect, c->G->A[node]);
	SetDelete(c->G->A[node], node);
	if(!SetEq(intersect, c->set))
	{
	    SetFree(intersect);
	    return NULL;
	}
    }
    assert(SetEq(intersect, c->set));
    SetFree(intersect);
    return c->set;
}

SET *GraphKnNext(CLIQUE *c)
{
    SET *s;
    while(CombinNext(c->combin))
	if((s = Graph_IsCombClique(c)))
	    return s;
    return NULL;
}

CLIQUE *GraphKnFirst(GRAPH *G, int k)
{
    int i, nDegk = G->n;
    SET *setDegk;
    CLIQUE *c;

    assert(k <= G->n);
    if(k == 0)
	return NULL;

    setDegk = SetAlloc(G->n);
    c = (CLIQUE*)Calloc(1,sizeof(CLIQUE));
    c->G = GraphCopy(NULL, G);
    c->cliqueSize = k;
    c->inducedArray = Calloc(G->n, sizeof(c->inducedArray[0]));

#if 1
    /*
    ** First reduce the potential number of vertices needed to check. A
    ** necessary condition for a vertex to be a member of a Kk is to have
    ** degree >= k-1, and all edges go to other vertices of degree >= k-1.
    ** Iterate inducing subgraphs 'til we cain't induce no more...
    */
    while(nDegk >= k)
    {
	int prevNDegk = nDegk;  /* just to check for changes */
	nDegk = 0;
	SetEmpty(setDegk);
	for(i=0; i < c->G->n; i++)
	    if(c->G->degree[i] >= k-1)
	    {
		++nDegk;
		SetAdd(setDegk, i);
	    }
	if(nDegk == prevNDegk)  /* nobody eliminated */
	    break;
	GraphInduced_NoVertexDelete(c->G, c->G, setDegk);
    }

    if(nDegk < k)
    {
	Free(c->inducedArray);
	GraphFree(c->G);
	Free(c);
	SetFree(setDegk);
	return NULL;
    }
    i = SetToArray(c->inducedArray, setDegk);
    assert(i == nDegk);
#else
    nDegk = G->n;
    for(i=0; i < nDegk; i++)
	c->inducedArray[i] = i;
#endif
    SetEmpty(setDegk);
    c->set = setDegk; /* re-use setDegK */
    c->combArray = Calloc(k, sizeof(c->combArray[0]));
    c->combin = CombinZeroth(nDegk, k, c->combArray);
    if(Graph_IsCombClique(c) || GraphKnNext(c))
	return c;
    /* else */
    GraphCliqueFree(c);
    return NULL;
}

void GraphCliqueFree(CLIQUE *c)
{
    CombinFree(c->combin);
    Free(c->combArray);
    Free(c->inducedArray);
    GraphFree(c->G);
    SetFree(c->set);
    Free(c);
}

CLIQUE *GraphInFirst(GRAPH *G, int n)
{
    static GRAPH Gbar;  /* non re-entrant */
    GraphComplement(&Gbar, G);
    return GraphKnFirst(&Gbar, n);
}


Boolean GraphKnContains(GRAPH *G, int n)
{
    CLIQUE *c;
    assert(n>=0);
    if(n < 2)
	return true;
    if(n > G->n)
	return false;
    c = GraphKnFirst(G, n);
    if(c)
	GraphCliqueFree(c);
    return c != NULL;
}

Boolean GraphInContains(GRAPH *G, int n)
{
    static GRAPH Gbar;  /* non re-entrant */
    GraphComplement(&Gbar, G);
    return GraphKnContains(&Gbar, n);
}

/**************************************************************************
**
**  Graph Isomorphism
**
**************************************************************************/

static GRAPH *isoG1, *isoG2;

static Boolean _permutationIdentical(int n, int perm[n])
{
    int i, j;
    for(i=0; i<n; i++)
	if(isoG1->degree[i] != isoG2->degree[perm[i]])
	    return false;

    for(i=0; i<n; i++) for(j=i+1; j<n; j++)
	/* The !GraphAreConnected is just to turn a bitstring into a boolean */
	if(!GraphAreConnected(isoG1, i,j) !=
	    !GraphAreConnected(isoG2, perm[i], perm[j]))
	    return false;   /* non-isomorphic */
    return true;   /* isomorphic! */
}

Boolean GraphsIsomorphic(int *perm, GRAPH *G1, GRAPH *G2)
{
    int i, n = G1->n, degreeCount1[n], degreeCount2[n];
    SET *degreeOnce;

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
    ** to the neighborhood of v2.
    */
    degreeOnce = SetAlloc(n);
    for(i=0; i<n; i++)
	if(degreeCount1[i] == 1)
	    SetAdd(degreeOnce, i);
    for(i=0; i<n; i++)
    {
	/* Find out if the degree of vertex i in G1 appears only once */
	if(SetIn(degreeOnce, G1->degree[i]))
	{
	    int j, degree = G1->degree[i];
	    GRAPH *neighG1i, *neighG2j;

	    /* find the (unique) vertex in G2 that has the same degree */
	    for(j=0; j < n; j++)
		if(G2->degree[j] == degree)
		    break;
	    assert(j < n);

	    assert((!G1->sparse||G1->sparse==both) && (!G2->sparse||G2->sparse==both));
	    neighG1i = GraphInduced(NULL, G1, G1->A[i]);
	    neighG2j = GraphInduced(NULL, G2, G2->A[j]);

	    /*
	    ** Note: this recursion works only as long as
	    ** _permutationIdentical doesn't call GraphsIsomorphic.
	    ** (if it does, isoG1 and isoG2 get messed up).
	    */

	    j = GraphsIsomorphic(perm, neighG1i, neighG2j);
	    GraphFree(neighG1i);
	    GraphFree(neighG2j);
	    if(!j)
		return false;
	    /* Otherwise they *might* be isomorphic, so keep going */
	}
    }
    SetFree(degreeOnce);

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
