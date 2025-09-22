#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "heap.h"
#include "queue.h"
#include "multisets.h"
#include "sorts.h"

#define VERBOSE 0
#define PARANOID_ASSERTS 1	// turn on paranoid checking --- slows down execution by a factor of 2-3
#define SPARSE true // do not try false at the moment, it's broken
static char *_BLANT_DIR=".";
#define CANON_DIR "canon_maps"

// Enable the code that uses C++ to parse input files?
static int *_pairs, _numNodes, _numEdges, _maxEdges=1024, _seed;
char **_nodeNames, _supportNodeNames = true;

#define USE_MarsenneTwister 0
#if USE_MarsenneTwister
#include "libwayne/MT19937/mt19937.h"
#define RandomSeed /*nothing*/
static MT19937 *_mt19937;
double RandomUniform(void) {
    if(!_mt19937) _mt19937 = Mt19937Alloc(_seed);
    return Mt19937NextDouble(_mt19937);
}
#else
#include "rand48.h"
#define RandomUniform drand48
#define RandomSeed srand48
#endif

int _sampleMethod = -1;
FILE *_sampleFile; // if _sampleMethod is SAMPLE_FROM_FILE
char *_sampleFileName; 
char _sampleFileEOF;

// Below are the sampling methods
#define SAMPLE_FROM_FILE 0
#define SAMPLE_ACCEPT_REJECT 1	// makes things REALLY REALLY slow.  Like 10-100 samples per second rather than a million.
#define SAMPLE_NODE_EXPANSION 2	// sample using uniform node expansion; about 100,000 samples per second
#define SAMPLE_EDGE_EXPANSION 3	// Fastest, up to a million samples per second
#define SAMPLE_RESERVOIR 4	// Lu Bressan's reservoir sampler, reasonably but not entirely unbiased.
#define SAMPLE_MCMC 5 // MCMC Algorithm estimates graphlet frequency with a random walk
#ifndef RESERVOIR_MULTIPLIER
// this*k is the number of steps in the Reservoir walk. 8 seems to work best, empirically.
#define RESERVOIR_MULTIPLIER 8
#endif

static Boolean _MCMC_UNIFORM = false; // Should MCMC restart at each edge

#define SAMPLE_FAYE 6

#define MAX_TRIES 100		// max # of tries in cumulative sampling before giving up

#define ALLOW_DISCONNECTED_GRAPHLETS 0

#define USE_INSERTION_SORT 0

// The following is the most compact way to store the permutation between a non-canonical and its canonical representative,
// when k=8: there are 8 entries, and each entry is a integer from 0 to 7, which requires 3 bits. 8*3=24 bits total.
// For simplicity we use the same 3 bits per entry, and assume 8 entries, even for k<8.  It wastes memory for k<4, but
// makes the coding much simpler.
typedef unsigned char kperm[3]; // The 24 bits are stored in 3 unsigned chars.

#define maxK 7
#define MAX_CANONICALS 1044
#define MAX_ORBITS 5096
#define maxBk (1U << (maxK*(maxK-1)/2))
#define LOWER_TRIANGLE 1
#define MAX_NODES 10000 // enough for HI-union from Vidal et al.

static unsigned int _Bk, _k; // _k is the global variable storing k; _Bk=actual number of entries in the canon_map for given k.
static unsigned _numCanon, _canonList[MAX_CANONICALS];
static SET *_connectedCanonicals;
static unsigned _numOrbits, _orbitList[MAX_CANONICALS][maxK];
static unsigned short _orbitCanonMapping[MAX_ORBITS]; // Maps orbits to canonical (including disconnected)

//windowRep global Variables
#define WINDOW_SAMPLE_MIN 1 // Find the k-graphlet with minimal canonicalInt
#define WINDOW_SAMPLE_MAX 2 // Find the k-graphlet with maximal canonicalInt
#define WINDOW_SAMPLE_MIN_D 3   // Find the k-graphlet with minimal canonicalInt and balanced numEdges
#define WINDOW_SAMPLE_MAX_D 4   // Find the k-graphlet with maximal canonicalInt and balanced numEdges
#define WINDOW_SAMPLE_LEAST_FREQ_MIN 5 // Find the k-graphlet with least fequent cacnonicalInt. IF there is a tie, pick the minimal one
#define WINDOW_SAMPLE_LEAST_FREQ_MAX 6 // Find the k-graphlet with least fequent cacnonicalInt. IF there is a tie, pick the maximial one 
int _windowSampleMethod = -1;

#define WINDOW_COMBO 0 // Turn on for using Combination method to sample k-graphlets in Window. Default is DFS-like way.
static int _windowSize = 0;
static Boolean _window = false;
static int** _windowReps;
static int _MAXnumWindowRep = 0;
static int _numWindowRep = 0;

enum OutputMode {undef, indexGraphlets, predictPairs, indexOrbits, graphletFrequency, outputODV, outputGDV};
static enum OutputMode _outputMode = undef;
static unsigned long int _graphletCount[MAX_CANONICALS];

enum CanonicalDisplayMode {undefined, ordinal, decimal, binary, orca, jesse};
static enum CanonicalDisplayMode _displayMode = undefined;

enum FrequencyDisplayMode {freq_display_mode_undef, count, concentration};
static enum FrequencyDisplayMode _freqDisplayMode = freq_display_mode_undef;

static int _magicTable[MAX_CANONICALS][12]; //Number of canonicals for k=8 by number of columns in magic table
static int _outputMapping[MAX_CANONICALS];

// A bit counter-intuitive: we need to allocate this many vectors each of length [_numNodes],
// and then the degree for node v, graphlet/orbit g is _degreeVector[g][v], NOT [v][g].
// We do this simply because we know the length of MAX_CANONICALS so we pre-know the length of
// the first dimension, otherwise we'd need to get more funky with the pointer allocation.
// Only one of these actually get allocated, depending upon outputMode.
static unsigned long int *_graphletDegreeVector[MAX_CANONICALS];
static unsigned long int    *_orbitDegreeVector[MAX_ORBITS];
static double *_doubleOrbitDegreeVector[MAX_ORBITS];

// If you're squeemish then use this one to access the degrees:
#define ODV(node,orbit)       _orbitDegreeVector[orbit][node]
#define GDV(node,graphlet) _graphletDegreeVector[graphlet][node]

// number of parallel threads to run.  This must be global because we may get called from C++.
static int _THREADS;

// Here's where we're lazy on saving memory, and we could do better.  We're going to allocate a static array
// that is big enough for the 256 million permutations from non-canonicals to canonicals for k=8, even if k<8.
// So we're allocating 256MBx3=768MB even if we need much less.  I figure anything less than 1GB isn't a big deal
// these days. It needs to be aligned to a page boundary since we're going to mmap the binary file into this array.
static kperm Permutations[maxBk] __attribute__ ((aligned (8192)));
// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
static unsigned short _K[maxBk] __attribute__ ((aligned (8192)));

//The number of edges required to walk a *Hamiltonion* path
static unsigned _MCMC_L; // walk length for MCMC algorithm. k-d+1 with d almost always being 2.
static int _alphaList[MAX_CANONICALS];
static double _graphletConcentration[MAX_CANONICALS];

/* AND NOW THE CODE */


// You provide a permutation array, we fill it with the permutation extracted from the compressed Permutation mapping.
// There is the inverse transformation, called "EncodePerm", in createBinData.c.
static void ExtractPerm(char perm[_k], int i)
{
    int j, i32 = 0;
    for(j=0;j<3;j++) i32 |= (Permutations[i][j] << j*8);
    for(j=0;j<_k;j++)
	perm[j] = (i32 >> 3*j) & 7;
}

/*
** Given an integer, build the graph into the TINY_GRAPH *G, which has already been allocated.
** Handles either upper or lower triangle representation depending upon compile-time option below.
*/
void Int2TinyGraph(TINY_GRAPH* G, unsigned Gint)
{
    int i, j, bitPos=0, k = G->n;
    unsigned Gint2 = Gint;  // Gint2 has bits nuked as they're used, so when it's zero we can stop.
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
	    unsigned bit = ((unsigned)1 << bitPos);
	    if(Gint & bit)
		TinyGraphConnect(G,i,j);
	    Gint2 &= ~bit;
	    bitPos++;
	    assert(bitPos < 8*sizeof(unsigned)); // technically they could be equal... change when that happens
	}
	if(!Gint2) break;
    }
}


// Given a TINY_GRAPH and k, return the integer ID created from one triangle (upper or lower) of the adjacency matrix.
unsigned TinyGraph2Int(TINY_GRAPH *g, int k)
{
    int i, j, bitPos=0;
    unsigned Gint = 0, bit;

#if LOWER_TRIANGLE	// Prefer lower triangle to be compatible with Ine Melckenbeeck's Jesse code.
    for(i=k-1;i>0;i--)
    {
        for(j=i-1;j>=0;j--)
#else   // UPPER_TRIANGLE // this is what we used in the original faye code and paper with Adib Hasan and Po-Chien Chung.
    for(i=k-2;i>=0;i--)
    {
        for(j=k-1;j>i;j--)
#endif
        {
	    if(TinyGraphAreConnected(g,i,j))
	    {
		bit = (((unsigned)1) << bitPos);
		Gint |= bit;
	    }
            bitPos++;
	    assert(bitPos < 8*sizeof(unsigned)); // technically they could be equal... change when that happens
        }
    }
    return Gint;
}

// Given the big graph G and a set of nodes in V, return the TINY_GRAPH created from the induced subgraph of V on G.
static TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, int *Varray)
{
    unsigned i, j;
    TinyGraphEdgesAllDelete(Gv);
    for(i=0; i < Gv->n; i++) for(j=i+1; j < Gv->n; j++)
        if(GraphAreConnected(G, Varray[i], Varray[j]))
            TinyGraphConnect(Gv, i, j);
    return Gv;
}


// return how many nodes found. If you call it with startingNode == 0 then we automatically clear the visited array
static int NumReachableNodes(TINY_GRAPH *g, int startingNode)
{
    static TSET _visited;
    if(startingNode == 0) TSetEmpty(_visited);
    TSetAdd(_visited,startingNode);
    int j, Varray[maxK], numVisited = 0;
    int numNeighbors = TSetToArray(Varray, g->A[startingNode]);
    assert(numNeighbors == g->degree[startingNode]);
    for(j=0; j<numNeighbors; j++)if(!TSetIn(_visited,Varray[j])) numVisited += NumReachableNodes(g,Varray[j]);
    return 1+numVisited;
}

static int _numConnectedComponents;
static int *_whichComponent; // will be an array of size G->n specifying which CC each node is in.
static int *_componentSize; // number of nodes in each CC
static int **_componentList; // list of lists of components, largest to smallest.
static double _totalCombinations, *_combinations, *_probOfComponent, *_cumulativeProb;
static SET **_componentSet;

void PrintNode(int v) {
    if(_supportNodeNames)
	printf("%s", _nodeNames[v]);
    else
	printf("%d", v);
}

static int InitializeConnectedComponents(GRAPH *G)
{
    static unsigned int v, *Varray, j, i;
    assert(!Varray); // we only can be called once.
    assert(_numConnectedComponents == 0);
    SET *visited = SetAlloc(G->n);
    Varray = Calloc(G->n, sizeof(int));
    _whichComponent = Calloc(G->n, sizeof(int));
    _componentSize = Calloc(G->n, sizeof(int)); // probably bigger than it needs to be but...
    _componentList = Calloc(G->n, sizeof(int*)); // probably bigger...
    _combinations = Calloc(G->n, sizeof(double*)); // probably bigger...
    _probOfComponent = Calloc(G->n, sizeof(double*)); // probably bigger...
    _cumulativeProb = Calloc(G->n, sizeof(double*)); // probably bigger...
    _componentSet = Calloc(G->n, sizeof(SET*));
    
    int nextStart = 0;
    _componentList[0] = Varray;
    for(v=0; v < G->n; v++) if(!SetIn(visited, v))
    {
	_componentSet[_numConnectedComponents] = SetAlloc(G->n);
	_componentList[_numConnectedComponents] = Varray + nextStart;
	GraphVisitCC(G, v, _componentSet[_numConnectedComponents], Varray + nextStart, _componentSize + _numConnectedComponents);
	SetUnion(visited, visited, _componentSet[_numConnectedComponents]);
	for(j=0; j < _componentSize[_numConnectedComponents]; j++)
	{
	    assert(_whichComponent[Varray[nextStart + j]] == 0);
	    _whichComponent[Varray[nextStart + j]] = _numConnectedComponents;
	}
	nextStart += _componentSize[_numConnectedComponents];
	++_numConnectedComponents;
    }
    assert(nextStart == G->n);

    _totalCombinations = 0.0;
    for(i=0; i< _numConnectedComponents; i++)
    {
	//find the biggest one
	int biggest = i;
	for(j=i+1; j<_numConnectedComponents;j++)
	    if(_componentSize[j] > _componentSize[biggest])
		biggest = j;
	// Now swap the biggest one into position i;
	for(j=0; j < _componentSize[biggest]; j++)
	    _whichComponent[_componentList[biggest][j]] = i;
	int itmp, *pitmp;
	SET * stmp;
	itmp = _componentSize[i];
	_componentSize[i] = _componentSize[biggest];
	_componentSize[biggest] = itmp;
	pitmp = _componentList[i];
	_componentList[i] = _componentList[biggest];
	_componentList[biggest] = pitmp;
    stmp = _componentSet[i];
    _componentSet[i] = _componentSet[biggest];
    _componentSet[biggest] = stmp;
	_combinations[i] = CombinChooseDouble(_componentSize[i], _k);
	_totalCombinations += _combinations[i];
    }

    double cumulativeProb = 0.0;
    for(i=0; i< _numConnectedComponents; i++)
    {
	_probOfComponent[i] =  _combinations[i] / _totalCombinations;
	_cumulativeProb[i] = cumulativeProb + _probOfComponent[i];
	cumulativeProb = _cumulativeProb[i];
	if(cumulativeProb > 1)
	{
	    assert(cumulativeProb - 1.0 < _numConnectedComponents * 1e-15);  // allow some roundoff error
	    cumulativeProb = _cumulativeProb[i] = 1.0;
	}
	//printf("Component %d has %d nodes and probability %lf, cumulative prob %lf\n", i, _componentSize[i], _probOfComponent[i], _cumulativeProb[i]);
    }
}

// Given the big graph G and an integer k, return a k-graphlet from G
// in the form of a SET of nodes called V. When complete, |V| = k.
// Caller is responsible for allocating the set V and its array Varray.
// The algorithm works by maintaining the exact set of edges that are
// emanating out of the graphlet-in-construction. Then we pick one of
// these edges uniformly at random, add the node at the endpoint of that
// edge to the graphlet, and then add all the outbound edges from that
// new node (ie., edges that are not going back inside the graphlet).
// So the "outset" is the set of edges going to nodes exactly distance
// one from the set V, as V is being built.

static SET *SampleGraphletNodeBasedExpansion(SET *V, int *Varray, GRAPH *G, int k, int whichCC)
{
    static SET *outSet;
    static int numIsolatedNodes;
    if(!outSet)
       outSet = SetAlloc(G->n);  // we won't bother to free this since it's static.
    else if(G->n > outSet->maxElem)
	SetResize(outSet, G->n);
    else
	SetEmpty(outSet);
    int v1, v2, i;
    int nOut = 0, outbound[G->n]; // vertices one step outside the boundary of V
    assert(V && V->maxElem >= G->n);
    SetEmpty(V);
    int edge;
    do {
	edge = G->numEdges * RandomUniform();
	v1 = G->edgeList[2*edge];
    } while(!SetIn(_componentSet[whichCC], v1));
    v2 = G->edgeList[2*edge+1];
    SetAdd(V, v1); Varray[0] = v1;
    SetAdd(V, v2); Varray[1] = v2;

    // The below loops over neighbors can take a long time for large graphs with high mean degree. May be faster
    // with bit operations if we stored the adjacency matrix... which may be too big to store for big graphs. :-(
    for(i=0; i < G->degree[v1]; i++)
    {
	int nv1 =  G->neighbor[v1][i];
	if(nv1 != v2)
	{
	    assert(!SetIn(V, nv1)); // assertion to ensure we're in line with faye
	    SetAdd(outSet, (outbound[nOut++] = nv1));
	}
    }
    for(i=0; i < G->degree[v2]; i++)
    {
	int nv2 =  G->neighbor[v2][i];
	if(nv2 != v1 && !SetIn(outSet, nv2))
	{
	    assert(!SetIn(V, nv2)); // assertion to ensure we're in line with faye
	    SetAdd(outSet, (outbound[nOut++] = nv2));
	}
    }
    for(i=2; i<k; i++)
    {
	int j;
	if(nOut == 0) // the graphlet has saturated it's connected component
	{
	    assert(SetCardinality(outSet) == 0);
	    assert(SetCardinality(V) < k);
#if ALLOW_DISCONNECTED_GRAPHLETS
	    while(SetIn(V, (j = G->n*RandomUniform())))
		; // must terminate since k <= G->n
	    outbound[nOut++] = j;
	    j = 0;
#else
	    static int depth;
	    depth++;
	    // must terminate eventually as long as there's at least one connected component with >=k nodes.
	    assert(depth < MAX_TRIES); // graph is too disconnected
	    V = SampleGraphletNodeBasedExpansion(V, Varray, G, k, whichCC);
	    depth--;
	    // Ensure the damn thing really *is* connected.
	    TINY_GRAPH *T = TinyGraphAlloc(k);
	    TinyGraphInducedFromGraph(T, G, Varray);
	    assert(NumReachableNodes(T,0) == k);
	    TinyGraphFree(T);
	    return V;
#endif
	}
	else
	    j = nOut * RandomUniform();
	v1 = outbound[j];
	SetDelete(outSet, v1);
	SetAdd(V, v1); Varray[i] = v1;
	outbound[j] = outbound[--nOut];	// nuke v1 from the list of outbound by moving the last one to its place
	for(j=0; j<G->degree[v1];j++) // another loop over neighbors that may take a long time...
	{
	    v2 = G->neighbor[v1][j];
	    if(!SetIn(outSet, v2) && !SetIn(V, v2))
		SetAdd(outSet, (outbound[nOut++] = v2));
	}
    }
    assert(i==k);
#if PARANOID_ASSERTS
    assert(SetCardinality(V) == k);
    assert(nOut == SetCardinality(outSet));
#endif
    return V;
}


// modelled after faye by Tuong Do
static SET *SampleGraphletFaye(SET *V, int *Varray, GRAPH *G, int k, int whichCC)
{
    /* Faye: Add a visited array to keep track of nodes. Initialize to 0 */
    int visited[G->n];
    static SET *outSet;
    static int numIsolatedNodes;
    if(!outSet)
       outSet = SetAlloc(G->n);  // we won't bother to free this since it's static.
    else if(G->n > outSet->maxElem)
	SetResize(outSet, G->n);
    else
	SetEmpty(outSet);
    int v1, v2, i;
    int nOut = 0, outbound[G->n]; // vertices one step outside the boundary of V
    assert(V && V->maxElem >= G->n);
    SetEmpty(V);
    int edge;
    do {
	edge = G->numEdges * RandomUniform();
	v1 = G->edgeList[2*edge];
    } while(!SetIn(_componentSet[whichCC], v1));
    v2 = G->edgeList[2*edge+1];
    SetAdd(V, v1); Varray[0] = v1;
    SetAdd(V, v2); Varray[1] = v2;

    /* Faye: Mark v1 and v2 as visited */
    visited[v1] = 1;
    visited[v2] = 1;    

    // The below loops over neighbors can take a long time for large graphs with high mean degree. May be faster
    // with bit operations if we stored the adjacency matrix... which may be too big to store for big graphs. :-(
    for(i=0; i < G->degree[v1]; i++)
    {
	int nv1 =  G->neighbor[v1][i];
	if(nv1 != v2)
	{
	    assert(!SetIn(V, nv1)); // assertion to ensure we're in line with faye
	    if (!visited[nv1]) { /* Faye: Check if it's visited */
            SetAdd(outSet, (outbound[nOut++] = nv1));
            visited[nv1] = 1;
        }
	}
    }
    for(i=0; i < G->degree[v2]; i++)
    {
	int nv2 =  G->neighbor[v2][i];
	if(nv2 != v1 && !SetIn(outSet, nv2))
	{
	    assert(!SetIn(V, nv2)); // assertion to ensure we're in line with faye
	    if (!visited[nv2]) { /* Faye: Check if it's visited */
            SetAdd(outSet, (outbound[nOut++] = nv2));
            visited[nv2] = 1;
        }
	}
    }
    for(i=2; i<k; i++)
    {
	int j;
	if(nOut == 0) // the graphlet has saturated it's connected component
	{
	    assert(SetCardinality(outSet) == 0);
	    assert(SetCardinality(V) < k);
#if ALLOW_DISCONNECTED_GRAPHLETS
	    /* Faye: check if the random node is visited instead 
        *while(SetIn(V, (j = G->n*RandomUniform()))    
        */
        while(visited[(j = G->n*RandomUniform())])
		; // must terminate since k <= G->n
	    outbound[nOut++] = j;
	    j = 0;
#else
	    static int depth;
	    depth++;
	    // must terminate eventually as long as there's at least one connected component with >=k nodes.
	    assert(depth < MAX_TRIES); // graph is too disconnected
	    V = SampleGraphletNodeBasedExpansion(V, Varray, G, k, whichCC);
	    depth--;
	    // Ensure the damn thing really *is* connected.
	    TINY_GRAPH *T = TinyGraphAlloc(k);
	    TinyGraphInducedFromGraph(T, G, Varray);
	    assert(NumReachableNodes(T,0) == k);
	    TinyGraphFree(T);
	    return V;
#endif
	}
	else
	    j = nOut * RandomUniform();
	v1 = outbound[j];
	SetDelete(outSet, v1);
	SetAdd(V, v1); Varray[i] = v1;
	outbound[j] = outbound[--nOut];	// nuke v1 from the list of outbound by moving the last one to its place
	for(j=0; j<G->degree[v1];j++) // another loop over neighbors that may take a long time...
	{
	    v2 = G->neighbor[v1][j];
        /* Faye: check if it's invisted instead
        * if(!SetIn(outSet, v2) && !SetIn(V, v2)) */
        if (!visited[v2]) {
		    SetAdd(outSet, (outbound[nOut++] = v2));
	        visited[v2] = 1;
        }
    }
    }
    assert(i==k);
#if PARANOID_ASSERTS
    assert(SetCardinality(V) == k);
    assert(nOut == SetCardinality(outSet));
#endif
    return V;
}

// Returns NULL if there are no more samples
static SET *SampleGraphletFromFile(SET *V, int *Varray, GRAPH *G, int k)
{
    SetEmpty(V);
	int i, numRead;
	char line[BUFSIZ];
	char *s = fgets(line, sizeof(line), _sampleFile);
	if(!s){
		_sampleFileEOF = 1; // forces exit below
		return NULL;
	}
	switch(k)
	{
	case 3: numRead = sscanf(line, "%d%d%d",Varray,Varray+1,Varray+2); break;
	case 4: numRead = sscanf(line, "%d%d%d%d",Varray,Varray+1,Varray+2,Varray+3); break;
	case 5: numRead = sscanf(line, "%d%d%d%d%d",Varray,Varray+1,Varray+2,Varray+3,Varray+4); break;
	case 6: numRead = sscanf(line, "%d%d%d%d%d%d",Varray,Varray+1,Varray+2,Varray+3,Varray+4,Varray+5); break;
	case 7: numRead = sscanf(line, "%d%d%d%d%d%d%d",Varray,Varray+1,Varray+2,Varray+3,Varray+4,Varray+5,Varray+6); break;
	case 8: numRead = sscanf(line, "%d%d%d%d%d%d%d%d",Varray,Varray+1,Varray+2,Varray+3,Varray+4,Varray+5,Varray+6,Varray+7); break;
	default: Fatal("unknown k value %d",k);
	}
	assert(numRead == k);
	for(k=0;i<k;i++){
		assert(Varray[i] >= 0 && Varray[i] < G->n);
		SetAdd(V, Varray[i]);
	}
    return V;
}


/*
** This is a faster graphlet sampling routine, although it may produce a
** distribution of graphlets that is further from "ideal" than the above.
** The difference is that in the above method we explicitly build and maintain
** the "outset" (the set of nodes one step outside V), but this can expensive
** when the mean degree of the graph G is high, since we have to add all the
** new edges emanating from each new node that we add to V.  Instead, here we
** are not going to explicitly maintain this outset.  Instead, we're simply
** going to remember the total number of edges emanating from every node in
** V *including* edges heading back into V, which are techically useless to us.
** This sum is just the sum of all the degrees of all the nodes in V.  Call this
** number "outDegree".  Then, among all the edges emanating out of every node in
** V, we pick a random node from V where the probablity of picking any node v is
** proportional to its degree.  Then we pick one of the edges emanating out of v
** uniformly at random.  Thus, every edge leaving every node in V has equal
** probability of being chosen---even though some of them lead back into V. If
** that is the case, then we throw away that edge and try the whole process again.
** The advantage here is that for graphs G with high mean degree M, the probability
** of picking an edge leading back into V is bounded by (k-1)/M, which becomes small
** as M gets large, and so we're very unlikely to "waste" samples.  Thus, unlike the
** above algorithm that gets *more* expensive with increasing density of G, this
** algorithm actually get cheaper with increasing density.
** Another way of saying this is that we avoid actually building the outSet
** as we do above, we instead simply keep a *estimated* count C[v] of separate
** outSets of each node v in the accumulating graphlet. Then we pick
** an integer in the range [0,..sum(C)), go find out what node that is,
** and if it turns out to be a node already in V, then we simply try again.
** It turns out that even if the mean degree is small, this method is *still*
** faster.  In other words, it's *always* faster than the above method.  The
** only reason we may prefer the above method is because it's theoretically cleaner
** to describe the distribution of graphlets that comes from it---although
** empirically this one does reasonably well too.  However, if the only goal
** is blinding speed at graphlet sampling, eg for building a graphlet database
** index, then this is the preferred method.
*/
static SET *SampleGraphletEdgeBasedExpansion(SET *V, int *Varray, GRAPH *G, int k, int whichCC)
{
    int edge, v1, v2;
    assert(V && V->maxElem >= G->n);
    SetEmpty(V);
    int nOut = 0;
    do {
	edge = G->numEdges * RandomUniform();
	v1 = G->edgeList[2*edge];
    } while(!SetIn(_componentSet[whichCC], v1));
    v2 = G->edgeList[2*edge+1];
    SetAdd(V, v1); Varray[0] = v1;
    SetAdd(V, v2); Varray[1] = v2;
    int vCount = 2;

    int outDegree = G->degree[v1] + G->degree[v2];
    static int cumulative[maxK];
    cumulative[0] = G->degree[v1]; // where v1 = Varray[0]
    cumulative[1] = G->degree[v2] + cumulative[0];

    static SET *internal;	// mark choices of whichNeigh that are discovered to be internal
    static int Gn;
    if(!internal) {internal = SetAlloc(G->n); Gn = G->n;}
    else if(Gn != G->n) {SetFree(internal); internal = SetAlloc(G->n); Gn=G->n;}
    else SetEmpty(internal);

    int numTries = 0;
    while(vCount < k)
    {
	int i, whichNeigh, newNode = -1;
	while(numTries < MAX_TRIES && SetIn(internal, (whichNeigh = outDegree * RandomUniform())))
	    ++numTries; // which edge to choose among all edges leaving all nodes in V so far?
	if(numTries >= MAX_TRIES) {
#if ALLOW_DISCONNECTED_GRAPHLETS
	    // get a new node outside this connected component.
	    // Note this will return a disconnected graphlet.
	    while(SetIn(V, (newNode = G->n*RandomUniform())))
		; // must terminate since k <= G->n
	    numTries = 0;
	    outDegree = 0;
	    int j;
	    for(j=0; j<vCount; j++)	// avoid picking these nodes ever again.
		cumulative[j] = 0;
	    SetEmpty(internal);
#else
	    static int depth;
	    depth++;
	    assert(depth < MAX_TRIES);
	    V = SampleGraphletEdgeBasedExpansion(V, Varray, G, k, whichCC);
	    depth--;
	    return V;
#endif
	}
	for(i=0; cumulative[i] <= whichNeigh; i++)
	    ; // figure out whose neighbor it is
	int localNeigh = whichNeigh-(cumulative[i]-G->degree[Varray[i]]); // which neighbor of node i?
	if(newNode < 0) newNode = G->neighbor[Varray[i]][localNeigh];
#if PARANOID_ASSERTS
	// really should check some of these a few lines higher but let's group all the paranoia in one place.
	assert(i < vCount);
	assert(0 <= localNeigh && localNeigh < G->degree[Varray[i]]);
	assert(0 <= newNode && newNode < G->n);
#endif
	if(SetIn(V, newNode))
	{
	    SetAdd(internal, whichNeigh);
	    if(++numTries < MAX_TRIES)
		continue;
	    else // graph is too disconnected
	    {
#if PARANOID_ASSERTS
		// We are probably in a connected component with fewer than k nodes.
		// Test that hypothesis.
		assert(!GraphCCatLeastK(G, v1, k));
#endif
#if ALLOW_DISCONNECTED_GRAPHLETS
		// get a new node outside this connected component.
		// Note this will return a disconnected graphlet.
		while(SetIn(V, (newNode = G->n*RandomUniform())))
		    ; // must terminate since k <= G->n
		numTries = 0;
		outDegree = 0;
		int j;
		for(j=0; j<vCount; j++)	// avoid picking these nodes ever again.
		    cumulative[j] = 0;
		SetEmpty(internal);
#else
		static int depth;
		depth++;
		assert(depth < MAX_TRIES);
		V = SampleGraphletEdgeBasedExpansion(V, Varray, G, k, whichCC);
		depth--;
		return V;
#endif
	    }
	}
	SetAdd(V, newNode);
	cumulative[vCount] = cumulative[vCount-1] + G->degree[newNode];
	Varray[vCount++] = newNode;
	outDegree += G->degree[newNode];
#if PARANOID_ASSERTS
	assert(SetCardinality(V) == vCount);
	assert(outDegree == cumulative[vCount-1]);
#endif
    }
#if PARANOID_ASSERTS
    assert(SetCardinality(V) == k);
    assert(vCount == k);
#endif
    return V;
}

// From the paper: ``Sampling Connected Induced Subgraphs Uniformly at Random''
// Xuesong Lu and Stephane Bressan, School of Computing, National University of Singapore
// {xuesong,steph}@nus.edu.sg
// 24th International Conference on Scientific and Statistical Database Management, 2012
// Note that they suggest *edge* based expansion to select the first k nodes, and then
// use reservoir sampling for the rest. But we know edge-based expansion sucks, so we'll
// start with a better starting guess, which is node-based expansion.
static SET *SampleGraphletLuBressanReservoir(SET *V, int *Varray, GRAPH *G, int k, int whichCC)
{
    // Start by getting the first k nodes using a previous method. Once you figure out which is
    // better, it's probably best to share variables so you don't have to recompute the outset here.
#if 1  // the following is copied almost verbatim from NodeEdgeExpansion, just changing for loop to while loop.
    static SET *outSet;
    static int numIsolatedNodes;
    if(!outSet)
       outSet = SetAlloc(G->n);  // we won't bother to free this since it's static.
    else if(G->n > outSet->maxElem)
	SetResize(outSet, G->n);
    else
	SetEmpty(outSet);
    int v1, v2, i;
    int nOut = 0, outbound[G->n]; // vertices one step outside the boundary of V
    assert(V && V->maxElem >= G->n);
    SetEmpty(V);
    int edge;
    do {
	edge = G->numEdges * RandomUniform();
	v1 = G->edgeList[2*edge];
    } while(!SetIn(_componentSet[whichCC], v1));
    v2 = G->edgeList[2*edge+1];
    SetAdd(V, v1); Varray[0] = v1;
    SetAdd(V, v2); Varray[1] = v2;

    // The below loops over neighbors can take a long time for large graphs with high mean degree. May be faster
    // with bit operations if we stored the adjacency matrix... which may be too big to store for big graphs. :-(
    for(i=0; i < G->degree[v1]; i++)
    {
	int nv1 =  G->neighbor[v1][i];
	if(nv1 != v2) SetAdd(outSet, (outbound[nOut++] = nv1));
    }
    for(i=0; i < G->degree[v2]; i++)
    {
	int nv2 =  G->neighbor[v2][i];
	if(nv2 != v1 && !SetIn(outSet, nv2)) SetAdd(outSet, (outbound[nOut++] = nv2));
    }
    i=2;

    // while(i<k || nOut > 0) // always do the loop at least k times, but i>=k is the reservoir phase.
    while(i < RESERVOIR_MULTIPLIER*k) // value of 8 seems to work best from empirical studies.
    {
	int candidate;
	if(nOut ==0) // the graphlet has saturated its connected component before getting to k, start elsewhere
	{
#if ALLOW_DISCONNECTED_GRAPHLETS
	    if(i < k)
	    {
		int tries=0;
		while(SetIn(V, (v1 = G->n*RandomUniform())))
		    assert(tries++<MAX_TRIES); // graph is too disconnected
		outbound[nOut++] = v1; // recall that nOut was 0 to enter this block, so now it's 1
		candidate = 0; // representing v1 as the 0'th entry in the outbound array
	    }
	    else
		assert(i==k); // we're done because i >= k and nOut == 0... but we shouldn't get here.
#else
	    static int depth;
	    depth++;
	    assert(depth < MAX_TRIES); // graph is too disconnected
	    V = SampleGraphletLuBressanReservoir(V, Varray, G, k, whichCC);
	    depth--;
	    return V;
#endif
	}
	else
	{
	    candidate = nOut * RandomUniform();
	    v1 = outbound[candidate];
	}
	assert(v1 == outbound[candidate]);
	SetDelete(outSet, v1);
	outbound[candidate] = outbound[--nOut];// nuke v1 by moving the last one to its place
	if(i < k)
	{
	    Varray[i] = v1;
	    SetAdd(V, v1);
	    int j;
	    for(j=0; j<G->degree[v1];j++) // another loop over neighbors that may take a long time...
	    {
		v2 = G->neighbor[v1][j];
		if(!SetIn(outSet, v2) && !SetIn(V, v2))
		    SetAdd(outSet, (outbound[nOut++] = v2));
	    }
	}
	else
	{
	    double reservoir_alpha = RandomUniform();
	    if(reservoir_alpha < k/(double)i)
	    {
		static TINY_GRAPH *T;
		if(!T) T = TinyGraphAlloc(k);
		static int graphetteArray[maxK], distArray[maxK];
#if PARANOID_ASSERTS
		// ensure it's connected before we do the replacement
		TinyGraphEdgesAllDelete(T);
		TinyGraphInducedFromGraph(T, G, Varray);
#if 0
		printf("NRN = %d\n", NumReachableNodes(T, 0));
		printf("BFS = %d\n", TinyGraphBFS(T, 0, k, graphetteArray, distArray));
#endif
		assert(NumReachableNodes(T, 0) == TinyGraphBFS(T, 0, k, graphetteArray, distArray));
		assert(NumReachableNodes(T, 0) == k);
#endif
		int memberToDelete = k*RandomUniform();
		v2 = Varray[memberToDelete]; // remember the node delated from V in case we need to revert
		Varray[memberToDelete] = v1; // v1 is the outbound candidate.
		TinyGraphEdgesAllDelete(T);
		TinyGraphInducedFromGraph(T, G, Varray);
		assert(NumReachableNodes(T,0) == TinyGraphBFS(T, 0, k, graphetteArray, distArray));
		if(NumReachableNodes(T, 0) < k)
		    Varray[memberToDelete] = v2; // revert the change because the graph is not connected
		else // add the new guy and delete the old
		{
#if PARANOID_ASSERTS
		    assert(SetCardinality(V) == k);
#endif
		    SetDelete(V, v2);
		    SetAdd(V, v1);
#if PARANOID_ASSERTS
		    assert(SetCardinality(V) == k);
#endif
		    int j;
		    for(j=0; j<G->degree[v1];j++) // another loop over neighbors that may take a long time...
		    {
			v2 = G->neighbor[v1][j];
			if(!SetIn(outSet, v2) && !SetIn(V, v2))
			    SetAdd(outSet, (outbound[nOut++] = v2));
		    }
		}
	    }
	}
	i++;
    }
#else
    SampleGraphletEdgeBasedExpansion(V, Varray, G, k, whichCC);
#endif
}

int alphaListPopulate(char *BUF, int *alpha_list, int k) {
	sprintf(BUF, "%s/%s/alpha_list_mcmc%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    int numAlphas, i;
    if(1!=fscanf(fp_ord, "%d",&numAlphas)) Fatal("read numAlphas failed") ;
    assert(numAlphas == _numCanon);
    for(i=0; i<numAlphas; i++) if(1!=fscanf(fp_ord, "%d", &alpha_list[i])) Fatal("read of alpha[%d] failed", i);
    fclose(fp_ord);
    return numAlphas;
}

#define mcmc_d 2

// Update the most recent d-graphlet to a random neighbor of it
int *MCMCGetNeighbor(int *Xcurrent, GRAPH *G)
{
	if (mcmc_d == 2)
	{
		int oldu = Xcurrent[0];
		int oldv = Xcurrent[1];
		int numTries = 0;
		while (oldu == Xcurrent[0] && oldv == Xcurrent[1]) {
#if PARANOID_ASSERTS
			assert(++numTries < MAX_TRIES);
#endif
			double p = RandomUniform();
			// if 0 < p < 1, p < deg(u) + deg(v) then
			if (p < ((double)G->degree[Xcurrent[0]])/(G->degree[Xcurrent[0]] + G->degree[Xcurrent[1]])) {
				// select randomly from Neigh(u) and swap
				int neighbor = (int) (G->degree[Xcurrent[0]] * RandomUniform());
				Xcurrent[1] = G->neighbor[Xcurrent[0]][neighbor];
			}
			else {
				// select randomly from Neigh(v) and swap
				int neighbor = (int) (G->degree[Xcurrent[1]] * RandomUniform());
				Xcurrent[0] = G->neighbor[Xcurrent[1]][neighbor];
			}
		}
#if PARANOID_ASSERTS
		assert(Xcurrent[0] != Xcurrent[1]);
		assert(oldu != Xcurrent[0] || oldv != Xcurrent[1]);
		assert(oldu != Xcurrent[1] || oldv != Xcurrent[0]);
#endif
	}
	else Fatal("Not implemented. Set d to 2");
	return Xcurrent;
}

// Crawls one step along the graph updating our sliding window
void crawlOneStep(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G) {
	int v, i;
	for (i = 0; i < mcmc_d; i++) { // Remove oldest d graphlet from sliding window
		v = QueueGet(XLQ).i;
		MultisetDelete(XLS, v);
	}
	MCMCGetNeighbor(X, G); // Gets a neighbor graphlet of the most recent d vertices and add to sliding window
	for (i = 0; i < mcmc_d; i++) {
		MultisetAdd(XLS, X[i]);
		QueuePut(XLQ, (foint) X[i]);
	}
}

// Initialize a sliding window represented by a multiset and queue of vertices.
// Sliding window is generated from a preselected edge from a predefined connected component and grown through edge walking.
void initializeSlidingWindow(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G, int windowSize, int edge)
{
	MultisetEmpty(XLS);
	QueueEmpty(XLQ);

	if (windowSize < 1) {
		Fatal("Window Size must be at least 1");
	}
#if PARANOID_ASSERTS
	assert(edge >= 0 && edge < G->numEdges);
#endif

	X[0] = G->edgeList[2 * edge];
	X[1] = G->edgeList[2 * edge + 1];

    MultisetAdd(XLS, X[0]); QueuePut(XLQ, (foint) X[0]);
    MultisetAdd(XLS, X[1]); QueuePut(XLQ, (foint) X[1]);
	
	// Add windowSize-1 d graphlets to our sliding window. The edge we added is the first d graphlet
	int i, j;
	for (i = 1; i < windowSize; i++) {
		MCMCGetNeighbor(X, G); // After each call latest graphlet is in X array
		for (j = 0; j < mcmc_d; j++) {
			MultisetAdd(XLS, X[j]);
			QueuePut(XLQ, (foint) X[j]);
		}
	}
}

// WalkLSteps fills XLS, XLQ (the sliding window) with L dgraphlets
// Given an empty sliding window, XLQ, XLS, walk along the graph starting at a random edge
// growing our sliding window until we have L graphlets in it.
// Then, we slide our window until it has k distinct vertices. That represents our initial sampling
// when we start/restart.
void WalkLSteps(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G, int k, int cc, int edge)
{
	
	//For now d must be equal to 2 because we start by picking a random edge
	if (mcmc_d != 2) {
		Fatal("mcmc_d must be set to 2 in blant.h for now");
	}

	if (edge < 0 && cc == -1) { // Pick a random edge from anywhere in the graph that has at least k nodes
		do {
		edge = G->numEdges * RandomUniform();
		X[0] = G->edgeList[2*edge];
		} while(!(_componentSize[_whichComponent[G->edgeList[2*edge]]] < k));
		X[1] = G->edgeList[2*edge+1];
	}
	else if (edge < 0) { // Pick a random edge from within a chosen connected component
		do {
		edge = G->numEdges * RandomUniform();
		X[0] = G->edgeList[2*edge];
		} while(!SetIn(_componentSet[cc], X[0]));
		X[1] = G->edgeList[2*edge+1];
	}
	// else start from the preselected edge

#if PARANOID_ASSERTS
	assert(_componentSize[_whichComponent[X[0]]] >= k); // Assert we can fill the window with the prechosen edge.
#endif

	initializeSlidingWindow(XLS, XLQ, X, G, _MCMC_L, edge);

	// Keep crawling until we have k distinct vertices
	static int numTries = 0;
	static int depth = 0;
	while (MultisetSupport(XLS) < k) {
		if (numTries++ > MAX_TRIES) { // If we crawl 100 steps without k distinct vertices restart
			assert(depth++ < MAX_TRIES); // If we restart 100 times in a row without success give up
			WalkLSteps(XLS,XLQ,X,G,k,cc,edge); // try again
			depth = 0; // If we get passed the recursive calls and successfully find a k graphlet, reset our depth
			numTries = 0; // And our number of attempts to crawl one step
			return;
		}
		crawlOneStep(XLS, XLQ, X, G);
	}
	numTries = 0;
}

static int _numSamples = 0;
static int _samplesPerEdge = 0;
/* SampleGraphletMCMC first starts with an edge in Walk L steps.
   It then walks along L = k-1 edges with MCMCGetNeighbor until it fills XLQ and XLS with L edges(their vertices are stored).
   XLQ and XLS always hold 2L vertices. They represent a window of the last L edges walked. If that window contains k
   distinct vertices, a graphlet is returned.
   This random walk predictably overcounts graphlets based on the alpha value and multipler.
   The alpha value is precomputed per graphlet type and the multiplier is based on the degree of the graphlets inside of it.
   After the algorithm is run the frequencies are normalized into concentrations.
*/

static SET *SampleGraphletMCMC(SET *V, int *Varray, GRAPH *G, int k, int whichCC) {
	static Boolean setup = false;
	static int currSamples = 0; // Counts how many samples weve done at the current starting point
	static int currEdge = 0; // Current edge we are starting at for uniform sampling
	static MULTISET *XLS = NULL; // A multiset holding L dgraphlets as separate vertex integers
	static QUEUE *XLQ = NULL; // A queue holding L dgraphlets as separate vertex integers
	static int Xcurrent[mcmc_d]; // holds the most recently walked d graphlet as an invariant
	static TINY_GRAPH *g = NULL; // Tinygraph for computing overcounting;
	if (!XLQ || !XLS || !g) {
		//NON REENTRANT CODE
		XLQ = QueueAlloc(k*mcmc_d);
		XLS = MultisetAlloc(G->n);
		g = TinyGraphAlloc(k);
	}

	// The first time we run this, or when we restart. We want to find our initial L d graphlets.
	if (!setup && !_MCMC_UNIFORM) {
		setup = true;
		WalkLSteps(XLS, XLQ, Xcurrent, G, k, whichCC, -1);
	}
	else if (_MCMC_UNIFORM && (!setup || currSamples >= _samplesPerEdge))
	{
		setup = true;
		WalkLSteps(XLS, XLQ, Xcurrent, G, k, whichCC, currEdge);
		do {
			currEdge++;
		} while (_componentSize[_whichComponent[G->edgeList[2*currEdge]]] < k);
		currSamples = 0;
	}
	else {
		// Keep crawling until we have k distinct vertices(can sample a graphlet). Crawl at least once
		do  {
			crawlOneStep(XLS, XLQ, Xcurrent, G);
		} while (MultisetSupport(XLS) != k);
		currSamples++;
	}
#if PARANOID_ASSERTS
		assert(MultisetSupport(XLS) == k); // very paranoid
		assert(QueueSize(XLQ) == 2 *_MCMC_L); // very paranoid
#endif

	/*
	Our sliding window now contains k distinct nodes. Fill the set V and array Varray with them
	The multiplier is a shorthand for d graphlet degree product, It is the product of the degrees
	of all the graphlets in the sliding window except the first and last.
	The degree of a graphlet is the sum of its outgoing edges.
	The alpha value is the number of ways to do a d-walk over the graphlet
	*/
	int node, numNodes = 0, i, j, graphletDegree;
	long long multiplier = 1;
	SetEmpty(V);

	for (i = 0; i < _MCMC_L; i++) {
		graphletDegree = -2; // The edge between the vertices in the graphlet isn't included and is double counted
		for (j = 0; j < mcmc_d; j++) {
			node = (XLQ->queue[(XLQ->front + (mcmc_d*i)+j) % XLQ->maxSize]).i;
			if (!SetIn(V, node)) {
				Varray[numNodes++] = node;
				SetAdd(V, node);
			}
			
			graphletDegree += G->degree[node];
		}
#if PARANOID_ASSERTS
		assert(graphletDegree > 0);
#endif

		// First and last graphlets in the window are skipped for multiplier product
		if (i != 0 && i != _MCMC_L-1) {
			multiplier *= (graphletDegree);
		}
		assert(multiplier > 0);
	}
	TinyGraphInducedFromGraph(g, G, Varray);
	int Gint = TinyGraph2Int(g, k);
	int GintCanon = _K[Gint];

#if PARANOID_ASSERTS
	assert(numNodes == k); // Ensure we are returning k nodes
#endif
	double count;
	if (_MCMC_L == 2) { // If _MCMC_L == 2, k = 3 and we can use the simplified overcounting formula.
		// The over counting ratio is the alpha value only.
		count += 1.0/(_alphaList[GintCanon]);
	} else {
		// The over counting ratio is the alpha value divided by the multiplier
		count += (double)multiplier/((double)_alphaList[GintCanon]);
	}
	if (_outputMode == outputODV) {
	    char perm[k];
	    memset(perm, 0, k);
	    ExtractPerm(perm, Gint);
	    for (j = 0; j < k; j++) {
		    _doubleOrbitDegreeVector[_orbitList[GintCanon][j]][Varray[(int)perm[j]]] += count;
	    }
	} else
	    _graphletConcentration[GintCanon] += count;

	return V; // return the sampled graphlet
}

// Fit SampleGraphletMCMC for windowRep implementation (alphalist and overcounting is not used here)
static SET* SampleWindowMCMC(SET *V, int *Varray, GRAPH *G, int W, int whichCC) 
{
	//Original SampleGraphletMCMC initial step. 
	// Not using tinyGraph to compute overcounting since W_size exceeeds the max tinygrpah size
	assert(W == _windowSize);
	static Boolean setup = false;
	static int currSamples = 0;
	static MULTISET *XLS = NULL; 
	static QUEUE *XLQ = NULL; 
	static int Xcurrent[mcmc_d]; 
	if (!XLQ || !XLS ) {
		XLQ = QueueAlloc(W*mcmc_d);
		XLS = MultisetAlloc(G->n);
	}

	if (!setup || (_numSamples/2 == currSamples++)) {
		setup = true;
		WalkLSteps(XLS, XLQ, Xcurrent, G, W, whichCC, -1);
	} else {
		do  {
			crawlOneStep(XLS, XLQ, Xcurrent, G);
		} while (MultisetSupport(XLS) != W);  //Keep crawling until we have W distinct vertices, Crawl at least once
	}
	int node, numNodes = 0, i, j, graphletDegree;
	SetEmpty(V);

	for (i = 0; i < _MCMC_L; i++) {
		graphletDegree = -2; //The edge between the vertices in the graphlet isn't included and is double counted
		for (j = 0; j < mcmc_d; j++) {
			node = (XLQ->queue[(XLQ->front + (mcmc_d*i)+j) % XLQ->maxSize]).i;
			if (!SetIn(V, node)) {
				Varray[numNodes++] = node;
				SetAdd(V, node);
			}
		}
	}
	return V;
}

// Convert the graphlet frequencies to concentrations
void finalizeMCMC() {
	double totalConcentration = 0;
	int i;
	for (i = 0; i < _numCanon; i++) {
		totalConcentration += _graphletConcentration[i];
#if PARANOID_ASSERTS
		assert(_graphletConcentration[i] >= 0.0);
#endif
	}
	for (i = 0; i < _numCanon; i++) {
		_graphletConcentration[i] /= totalConcentration;
	}
}

// Loads alpha values(overcounting ratios) for MCMC sampling from files
// The alpha value represents the number of ways to walk over that graphlet
// Global variable _MCMC_L is needed by many functions
// _MCMC_L represents the length of the sliding window in d graphlets for sampling
// Global variable _numSamples needed for the algorithm to reseed halfway through
// Concentrations are initialized to 0
void initializeMCMC(GRAPH* G, int k, int numSamples) {
	_MCMC_L = k - mcmc_d  + 1;
	// Count the number of valid edges to start from
	int i, validEdgeCount = 0;
	for (i = 0; i < G->numEdges; i++)
		if (_componentSize[_whichComponent[G->edgeList[2*i]]] >= k)
			validEdgeCount++;
	_samplesPerEdge =  (numSamples + (validEdgeCount / 2)) / validEdgeCount; // Division rounding up for samples per edge

	char BUF[BUFSIZ];
	_numSamples = numSamples;
	if(!_window)
	{
		alphaListPopulate(BUF, _alphaList, k);
		for (i = 0; i < _numCanon; i++)
		{
			_graphletConcentration[i] = 0.0;
		}
	}
}

// Compute the degree of the state in the state graph (see Lu&Bressen)
// Given the big graph G, and a set of nodes S (|S|==k), compute the 
// degree of the *state* represented by these nodes, which is:
//     Degree(S) = k*|NN2| + (k-1)*|NN1|
// where NN1 is the set of nodes one step outside S that have only 1 connection back into S,
// and   NN2 is the set of nodes one step outside S that have more than 1 connection back into S.
static int StateDegree(GRAPH *G, SET *S)
{
#if PARANOID_ASSERTS
    assert(SetCardinality(S) == _k);
#endif
    int Varray[_k]; // array of elements in S
    SetToArray(Varray, S);

    SET *outSet=SetAlloc(G->n); // the set of nodes in G that are one step outside S, from *anybody* in S
    int connectCount[G->n]; // for each node in the outset, count the number of times it's connected into S
    memset(connectCount, 0, G->n * sizeof(int)); // set them all to zero

    int i, j, numOut = 0;
    for(i=0;i<_k;i++) for(j=0; j<G->degree[Varray[i]]; j++)
    {
	int neighbor = G->neighbor[Varray[i]][j];
	if(!SetIn(S, neighbor))
	{
	    SetAdd(outSet, neighbor);
	    if(connectCount[neighbor] == 0) ++numOut;
	    ++connectCount[neighbor];
	}
    }
#if PARANOID_ASSERTS
    assert(SetCardinality(outSet) == numOut);
#endif
    int outArray[numOut], Sdeg = 0;
    i = SetToArray(outArray, outSet);
    assert(i == numOut);
    for(i=0; i < numOut; i++)
    {
	switch(connectCount[outArray[i]])
	{
	case 0: break;
	case 1: Sdeg += _k-1; break;
	default:  Sdeg += _k; break;
	}
    }
    SetFree(outSet);
    return Sdeg;
}


static SET *SampleGraphletLuBressan_MCMC_MHS_without_Ooze(SET *V, int *Varray, GRAPH *G, int k) { } // slower
static SET *SampleGraphletLuBressan_MCMC_MHS_with_Ooze(SET *V, int *Varray, GRAPH *G, int k) { } // faster!


/*
* Very slow: sample k nodes uniformly at random and throw away ones that are disconnected.
*/
static unsigned long int _acceptRejectTotalTries;
static SET *SampleGraphletAcceptReject(SET *V, int *Varray, GRAPH *G, int k)
{
    int arrayV[k], i;
    int nodeArray[G->n], distArray[G->n];
    TINY_GRAPH *g = TinyGraphAlloc(k);
    int graphetteArray[k];

    int tries = 0;
    do
    {
	SetEmpty(V);
	// select k nodes uniformly at random from G without regard to connectivity
	for(i=0; i<k; i++)
	{
	    do
		Varray[i] = G->n * RandomUniform();
	    while(SetIn(V, Varray[i]));
	    SetAdd(V, Varray[i]);
	}
	TinyGraphEdgesAllDelete(g);
	TinyGraphInducedFromGraph(g, G, Varray);
#if PARANOID_ASSERTS
	assert(SetCardinality(V)==k);
	assert(NumReachableNodes(g,0) == TinyGraphBFS(g, 0, k, graphetteArray, distArray));
#endif
	++tries;
    } while(NumReachableNodes(g, 0) < k);
    _acceptRejectTotalTries += tries;
    return V;
}


static int IntCmp(const void *a, const void *b)
{
    int *i = (int*)a, *j = (int*)b;
    return (*i)-(*j);
}

SET *canonListPopulate(char *BUF, int *canon_list, int k) {
    sprintf(BUF, "%s/%s/canon_list%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    unsigned short numCanon=0, i;
    int connected;
    if(1!=fscanf(fp_ord, "%hu\n",&numCanon) || numCanon==0) Fatal("canonListPopulate failed to read numCanon");
    SET *connectedCanonicals = SetAlloc(numCanon);
    for(i=0; i<numCanon; i++) {
	char buf[BUFSIZ], *tmp;
	tmp = fgets(buf, sizeof(buf), fp_ord);
	assert(tmp == buf);
	int len = strlen(buf);
	assert(buf[len] == '\0' && buf[len-1] == '\n');
	unsigned char canon_num_edges; // don't need this
	if(3!=sscanf(buf, "%u\t%d %hhd", &canon_list[i], &connected, &canon_num_edges))
	    Fatal("canonListPopulate: failed to read 3 inputs from line %d", i);
	if(connected) SetAdd(connectedCanonicals, i);
    }
    fclose(fp_ord);
    return connectedCanonicals;
}

unsigned orbitListPopulate(char *BUF,
	unsigned orbit_list[MAX_CANONICALS][maxK],
	unsigned short orbit_canon_mapping[MAX_ORBITS],
	unsigned short numCanon, int k) {
    sprintf(BUF, "%s/%s/orbit_map%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    unsigned o, numOrbits;
    if(1!=fscanf(fp_ord, "%u", &numOrbits)) Fatal("orbitListPopulate failed to read numOrbits");
    unsigned short c; int j;
    for(c=0; c<numCanon; c++) for(j=0; j<k; j++) {
	if(1!=fscanf(fp_ord, "%u", &orbit_list[c][j])) Fatal("orbitListPopulate: failed to read canon %d", c);
	orbit_canon_mapping[orbit_list[c][j]] = c;
    }
    fclose(fp_ord);
    return numOrbits;
}

/*
** Given a pre-allocated filename buffer, a 256MB aligned array K, num nodes k
** Mmap the canon_map binary file to the aligned array.
*/
unsigned short* mapCanonMap(char* BUF, unsigned short *K, int k) {
    int Bk = (1U <<(k*(k-1)/2));
    sprintf(BUF, "%s/%s/canon_map%d.bin", _BLANT_DIR, CANON_DIR, k);
    int Kfd = open(BUF, 0*O_RDONLY);
    if(Kfd <= 0) return NULL;
    unsigned short *Kf = Mmap(K, Bk*sizeof(short int), Kfd); // Using Mmap may cause error due to MAP_FIXED flag
    //unsigned short *Kf = (unsigned short*) mmap(K, sizeof(unsigned short)*Bk, PROT_READ, MAP_PRIVATE, Kfd, 0);
    assert(Kf != MAP_FAILED);
    return Kf;
}


// Assuming the global variable _k is set properly, go read in and/or mmap the big global
// arrays related to canonical mappings and permutations.
void SetGlobalCanonMaps(void)
{
    int i;
    char BUF[BUFSIZ];
    assert(3 <= _k && _k <= 8);
    _Bk = (1 <<(_k*(_k-1)/2));
    _connectedCanonicals = canonListPopulate(BUF, _canonList, _k);
    _numCanon = _connectedCanonicals->maxElem;
    _numOrbits = orbitListPopulate(BUF, _orbitList, _orbitCanonMapping, _numCanon, _k);
    void *v = mapCanonMap(BUF, _K, _k);
    //fprintf(stderr, "mapCanonMap %x, _K %x\n", v, _K);
    sprintf(BUF, "%s/%s/perm_map%d.bin", _BLANT_DIR, CANON_DIR, _k);
    int pfd = open(BUF, 0*O_RDONLY);
    kperm *Pf = Mmap(Permutations, _Bk*sizeof(Permutations[0]), pfd);
    assert(Pf == Permutations);
}

void LoadMagicTable()
{
	int i,j;
	char BUF[BUFSIZ];
	sprintf(BUF, "%s/orca_jesse_blant_table/UpperToLower%d.txt", _BLANT_DIR, _k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
	for(i=0; i<_numCanon; i++) {
	for(j=0; j<12 ;j++) {
	if(1!=fscanf(fp_ord, "%d", &_magicTable[i][j]))Fatal("magicTable read failed");
	}}
    fclose(fp_ord);
	switch (_displayMode) {
	case orca: // Load mapping from lower_ordinal to ORCA/Jesse ID into table for fast lookup
	case jesse:
	for (i=0; i < _numCanon; i++) {
		_outputMapping[_magicTable[i][4]] = _magicTable[i][7];
	}
	break;
	}
}

void SetBlantDir() {
    char* temp = getenv("BLANT_DIR");
    if (temp != NULL)
	_BLANT_DIR = strdup(temp); // can't assume the string returned by getetv never changes, so copy it.
}

void SampleGraphlet(GRAPH *G, SET *V, unsigned Varray[], int k) {
    int cc;
    double randomComponent = RandomUniform();
    for(cc=0; cc<_numConnectedComponents;cc++)
	if(_cumulativeProb[cc] > randomComponent)
	    break;
    //printf("choosing from CC %d\n", cc);
    switch (_sampleMethod) {
    case SAMPLE_ACCEPT_REJECT:
	SampleGraphletAcceptReject(V, Varray, G, k);	// REALLY REALLY SLOW and doesn't need to use cc
	break;
    case SAMPLE_NODE_EXPANSION:
	SampleGraphletNodeBasedExpansion(V, Varray, G, k, cc);
	break;
    case SAMPLE_FAYE:
	SampleGraphletFaye(V, Varray, G, k, cc);
	break;
    case SAMPLE_RESERVOIR:
	SampleGraphletLuBressanReservoir(V, Varray, G, k, cc); // pretty slow but not as bad as unbiased
	break;
    case SAMPLE_EDGE_EXPANSION:
	SampleGraphletEdgeBasedExpansion(V, Varray, G, k, cc); // Faster than NBE but less well tested and understood.
	break;
    case SAMPLE_MCMC:
	if(!_window)
		SampleGraphletMCMC(V, Varray, G, k, cc);
	else
		SampleWindowMCMC(V, Varray, G, k, cc);
	break;
    case SAMPLE_FROM_FILE:
	SampleGraphletFromFile(V, Varray, G, k);
	break;
    default:
	Fatal("unknown sampling method");
    }
}

void PrintCanonical(int GintCanon)
{
    int j, GintNumBits = _k*(_k-1)/2;
    char GintBinary[GintNumBits+1]; //Only used in -db output mode for indexing
    switch (_displayMode) {
    case undefined:
    case ordinal:
	printf("%d", GintCanon); // Note this is the ordinal of the canonical, not its bit representation
	break;
    case decimal: //Prints the decimal integer form of the canonical
	printf("%d", _canonList[GintCanon]);
	break;
    case binary: //Prints the bit representation of the canonical
	for (j=0;j<GintNumBits;j++)
	    {GintBinary[GintNumBits-j-1]=(((unsigned)_canonList[GintCanon] >> j) & 1 ? '1' : '0');}
	GintBinary[GintNumBits] = '\0';
	printf("%s", GintBinary);
	break;
	case orca: //Prints the ORCA ID of the canonical. Jesse uses same number.
	case jesse:
	printf("%d", _outputMapping[GintCanon]);
	break;
    }
}

void VarraySort(int *Varray, int k)
{
#if USE_INSERTION_SORT
    InsertionSortInt(Varray,k);
    //InsertionSort((void*)Varray,k,sizeof(Varray[0]),IntCmp);
#else
    qsort((void*)Varray, k, sizeof(Varray[0]), IntCmp);
#endif
}

static int _predictOrd = -1;
static int _predictOrbit1 = -1, _predictOrbit2=-1; // the columns (ie orbits) you want
static int _PredictMotifCount[MAX_CANONICALS][maxK][maxK];
static float _PredictCount[MAX_NODES][MAX_NODES];
static GRAPH *_G;

// Recursively count the specified orbits in the specified canonical motifs of g.
void PreProcessPredict(TINY_GRAPH *g, int topCanon)
{
    static int depth;
#if PARANOID_ASSERTS
    assert(g->n == _k && NumReachableNodes(g, 0)==g->n);
#endif
    int i,j, Gint = TinyGraph2Int(g,_k), GintOrd=_K[Gint];
    char perm[maxK];
    memset(perm, 0, _k);
    ExtractPerm(perm, Gint);
#if PARANOID_ASSERTS
    if(depth == 0) {
	assert(topCanon == GintOrd); // we should only be called on canonicals
	for(i=0;i<_k;i++) assert(perm[i]==i);
    }
#endif
    // Check to see if the whole thing is the cononical of interest before we start removing edges.
    if(GintOrd == _predictOrd) {
	for(i=0;i<_k-1;i++) for(j=i+1;j<_k;j++) if(!TinyGraphAreConnected(g,perm[i],perm[j])) {
	    if( (_orbitList[GintOrd][i] == _predictOrbit1 && _orbitList[GintOrd][j] == _predictOrbit2)
	     || (_orbitList[GintOrd][i] == _predictOrbit2 && _orbitList[GintOrd][j] == _predictOrbit1))
		++_PredictMotifCount[topCanon][(int)perm[i]][(int)perm[j]];
	}
    }

    // Now go about deleting edges recursively.
    for(i=0; i<_k-1; i++)for(j=i+1;j<_k;j++)
    {
	if(TinyGraphAreConnected(g,i,j)) // if it's an edge, delete it.
	{
	    TinyGraphDisconnect(g,i,j);
	    if(NumReachableNodes(g,0)==g->n) {
		++depth;
		PreProcessPredict(g,topCanon);
		--depth;
	    }
	    TinyGraphConnect(g,i,j);
	}
    }
}

void Predict_ProcessGraphlet(TINY_GRAPH *g, unsigned Varray[])
{
    static int depth;
#if PARANOID_ASSERTS
    assert(g->n == _k);
    assert(NumReachableNodes(g, 0)==g->n);
#endif
    int i,j, Gint = TinyGraph2Int(g,_k), GintOrd=_K[Gint];
    char perm[maxK];
    memset(perm, 0, _k);
    ExtractPerm(perm, Gint);
    for(i=0;i<_k-1;i++) for(j=i+1;j<_k;j++) {
	int g_u=perm[i], g_v=perm[j];
	int G_u=Varray[g_u], G_v=Varray[g_v];
	_PredictCount[G_u][G_v] += (_PredictMotifCount[GintOrd][i][j] + _PredictMotifCount[GintOrd][j][i]);
    }
}

void ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], char perm[], TINY_GRAPH *g, int k) {
    // We should probably figure out a faster sort? This requires a function call for every comparison.
    TinyGraphInducedFromGraph(g, G, Varray);
    int Gint = TinyGraph2Int(g,k), j, GintCanon=_K[Gint];
#if PARANOID_ASSERTS
    assert(0 <= GintCanon && GintCanon < _numCanon);
#endif
    switch(_outputMode)
    {
	static SET* printed;
    case graphletFrequency:
	++_graphletCount[GintCanon];
	break;
    case indexGraphlets:
#if SORT_INDEX_MODE // Note this destroys the columns-are-identical property, don't use by default.
	VarraySort(Varray, k);
#endif
	memset(perm, 0, k);
	ExtractPerm(perm, Gint);
	PrintCanonical(GintCanon);
	for(j=0;j<k;j++) {printf(" "); PrintNode(Varray[(int)perm[j]]);}
	puts("");
	break;
    case predictPairs:
	Predict_ProcessGraphlet(g,Varray);
	break;
    case indexOrbits:
#if SORT_INDEX_MODE // Note this destroys the columns-are-identical property, don't use by default.
	VarraySort(Varray, k);
#endif
	if(!printed) printed = SetAlloc(_k);
	SetEmpty(printed);
	memset(perm, 0, k);
	ExtractPerm(perm, Gint);
	PrintCanonical(GintCanon);
	for(j=0;j<k;j++) if(!SetIn(printed,j))
	{
	    putchar(' '); PrintNode(Varray[(int)perm[j]]);
	    SetAdd(printed, j);
	    int j1;
	    for(j1=j+1;j1<k;j1++) if(_orbitList[GintCanon][j1] == _orbitList[GintCanon][j])
	    {
		assert(!SetIn(printed, j1));
		putchar(':'); PrintNode(Varray[(int)perm[j1]]);
		SetAdd(printed, j1);
	    }
	}
	putchar('\n');
	break;
    case outputGDV:
	for(j=0;j<k;j++) ++GDV(Varray[j], GintCanon);
	break;
    case outputODV:
	memset(perm, 0, k);
	ExtractPerm(perm, Gint);
#if PERMS_CAN2NON            
	for(j=0;j<k;j++) ++ODV(Varray[(int)perm[j]], _orbitList[GintCanon][j]);
#else
	for(j=0;j<k;j++) ++ODV(Varray[j], _orbitList[GintCanon][(int)perm[j]]);
#endif
	break;
	    
    default: Abort("unknown or un-implemented outputMode");
	break;
    }
}


// Self construct the adjacency matrix of the window. Use _connectedCanonicals to check connectivity of the K-node graphlet
// No need to use TinyGraphInducedFromGraph, expensive calling GraphAreConnected for each combination
// This method is twice faster than previous
int combWindow2Int(int (*windowAdjList)[_windowSize], int *Varray, int *numEdges)
{
    int i, j, bitPos=0, Gint = 0, bit;
    *numEdges = 0;
    for(i=_k-1; i>0; i--)
        for(j=i-1;j>=0;j--)
        {
            if (windowAdjList[Varray[i]][Varray[j]] == 1)
            {
                bit = (1 << bitPos);
                Gint |= bit;
                *numEdges = *numEdges + 1;
            }
            bitPos++;
        }
    return Gint;
}

int getMaximumIntNumber(int K) 
{
    assert(K >= 3 && K <= 8);
    int num_of_bits = K * (K-1) / 2;
    return pow(2, num_of_bits);
}

int getD(int num_of_edges) 
{
    int M = _k * (_k - 1) / 2;
    int D = abs(2 * num_of_edges - M);
    return D;
}


void updateWindowRep(int *windowRepInt, int *D, int GintCanon, int pending_D, int *WArray, int *VArray, MULTISET *canonMSET)
{
    int i;
    if (_windowSampleMethod == WINDOW_SAMPLE_MIN || _windowSampleMethod == WINDOW_SAMPLE_MAX)
    {
        if(_windowSampleMethod == WINDOW_SAMPLE_MIN)
            if(GintCanon < *windowRepInt) {*windowRepInt = GintCanon; _numWindowRep = 0;}
        if(_windowSampleMethod == WINDOW_SAMPLE_MAX)
            if(GintCanon > *windowRepInt) {*windowRepInt = GintCanon; _numWindowRep = 0;}
        if(GintCanon == *windowRepInt)
        {
            for(i=0; i<_k; i++) _windowReps[_numWindowRep][i] = WArray[VArray[i]];
            _numWindowRep++;
        }
    }
    else if (_windowSampleMethod == WINDOW_SAMPLE_MIN_D || _windowSampleMethod == WINDOW_SAMPLE_MAX_D)
    {
        if (_windowSampleMethod == WINDOW_SAMPLE_MIN_D) 
            if(pending_D < *D || pending_D == *D && GintCanon < *windowRepInt) {*windowRepInt = GintCanon; *D = pending_D; _numWindowRep = 0;}
        if (_windowSampleMethod == WINDOW_SAMPLE_MAX_D)
            if(pending_D < *D || pending_D == *D && GintCanon > *windowRepInt) {*windowRepInt = GintCanon; *D = pending_D; _numWindowRep = 0;}
        if(pending_D == *D && GintCanon == *windowRepInt)
        {
            for(i=0; i<_k; i++) _windowReps[_numWindowRep][i] = WArray[VArray[i]];
            _numWindowRep++;
        }
    }
    else if (_windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MIN || _windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MAX)
    {
        for(i=0; i<_k; i++) _windowReps[_numWindowRep][i] = WArray[VArray[i]];
        _windowReps[_numWindowRep][_k] = GintCanon;
        if(canonMSET->array[GintCanon] < MAX_MULTISET_FREQ) MultisetAdd(canonMSET, GintCanon);
        _numWindowRep++;
    }
    else
        Fatal("unknown window sampling method.");
}

void updateLeastFrequent(int *windowRepInt, MULTISET *canonMSET) 
{
    int ordinals[MultisetSupport(canonMSET)];
    int i, pos = SetToArray(ordinals, canonMSET->set);
    int freq = _numWindowRep, multiplicity;
    for(i = 0; i < pos; i++)
    {
        multiplicity = MultisetMultiplicity(canonMSET, ordinals[i]);
        if(multiplicity == freq) {
            if (_windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MIN)
                if(ordinals[i] < *windowRepInt) *windowRepInt = ordinals[i];
            if (_windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MAX)
                if(ordinals[i] > *windowRepInt) *windowRepInt = ordinals[i]; 
        }
        else if(multiplicity < freq) {
            *windowRepInt = ordinals[i];
            freq = multiplicity;
        }
    }
    int new_numWindowRep = 0, j;
    for(i=0; i < _numWindowRep; i++)
    {
        if(_windowReps[i][_k] == *windowRepInt) {
            for(j=0; j<_k; j++) _windowReps[new_numWindowRep][j] = _windowReps[i][j];
            new_numWindowRep++;
        }
    }
    _numWindowRep = new_numWindowRep;
}

        
void ExtendSubGraph(GRAPH *Gi, int *WArray, int *VArray, SET *Vextension, int v, int *varraySize, int(*windowAdjList)[_windowSize], int *windowRepInt, int *D, MULTISET *canonMSET)
{
    int u, w, i, j, Gint, GintCanon, GintCanonInt, pending_D, numEdges=0;
    Boolean inclusive = false;
    SET *Vext = SetAlloc(Gi->n), *uNeighbors = SetAlloc(Gi->n);
    if(*varraySize == _k)
    {
        Gint = combWindow2Int(windowAdjList, VArray, &numEdges);
        GintCanon = _K[Gint];
        GintCanonInt = _canonList[GintCanon];
        pending_D = getD(numEdges);
        updateWindowRep(windowRepInt, D, GintCanon, pending_D, WArray, VArray, canonMSET);
    }
    else
    {
        while(SetCardinality(Vextension) != 0)
        {
            SetEmpty(Vext);
            // Remove an element w from Vextension
            w = Vextension->smallestElement;
            SetDelete(Vextension, (int)w);
            // Add exlusive neighbor u of w and u > v
            for(i=0; i<Gi->degree[w]; i++)
            {
                u = Gi->neighbor[w][i]; inclusive = false;
                if(u > v)
                {
                    SetEmpty(uNeighbors);
                    SetFromArray(uNeighbors, Gi->degree[u], Gi->neighbor[u]);
                    for(j=0; j<*varraySize; j++) if(SetIn(uNeighbors, VArray[j])) {inclusive = true; break;}
                    if(!inclusive && u > v) SetAdd(Vext, u);
                }
            }
            int* VArrayCopy = Calloc(_k, sizeof(int));
            for(i=0; i<_k; i++) VArrayCopy[i] = VArray[i];
            int varrayCopySize = *varraySize;
            VArrayCopy[varrayCopySize++] = w;
            SetUnion(Vext, Vext, Vextension);
            ExtendSubGraph(Gi, WArray, VArrayCopy, Vext, v, &varrayCopySize, windowAdjList, windowRepInt, D, canonMSET);
            free(VArrayCopy);
        }
    }
    SetFree(Vext);
    SetFree(uNeighbors);
    return;
}

// Right now use least frequent windowRep canonicals 
void FindWindowRepInWindow(GRAPH *G, SET *W, int *windowRepInt, int *D)
{
    int WArray[_windowSize], *VArray, ca[_k], i, j, Gint, GintCanon, GintCanonInt, numEdges=0, pending_D;   // Window node array, pending_window node array
    assert(SetToArray(WArray, W) == _windowSize);
    MULTISET *canonMSET = MultisetAlloc(getMaximumIntNumber(_k));
    COMBIN *c = CombinZeroth(_windowSize, _k, ca);  // (W choose K) many k-node graphlets in Window
    
    // Construct Adj Matrix for the Window
    int windowAdjList[_windowSize][_windowSize];
    for(i=0; i<_windowSize;i++)
        for(j=i+1;j<_windowSize;j++)
        {
            if(GraphAreConnected(G, WArray[i], WArray[j]))
                {windowAdjList[i][j]=1; windowAdjList[j][i]=1;}
            else
                {windowAdjList[i][j]=0; windowAdjList[j][i]=0;}
        }
        
    // Sampling K-graphlet Step
    VArray = Calloc(_k, sizeof(int));
    if(WINDOW_COMBO) 
    {
        do
        {
            for(i=0; i<_k; i++) 
                VArray[i] = ca[i];
            Gint = combWindow2Int(windowAdjList, VArray, &numEdges);
            GintCanon = _K[Gint];
            if(SetIn(_connectedCanonicals, GintCanon))
            {
                GintCanonInt = _canonList[GintCanon];
                pending_D = getD(numEdges);
                updateWindowRep(windowRepInt, D, GintCanon, pending_D, WArray, VArray, canonMSET);
            } 
        } while(CombinNext(c));
    }
    else 
    {
        GRAPH *Gi = GraphInduced(G, W);
        int v, varraySize;
        SET *Vextension = SetAlloc(Gi->n);
        for(v=0; v<_windowSize; v++)
        {
            SetEmpty(Vextension);
            varraySize=0;
            for(i=0; i<Gi->degree[v]; i++) 
                if(Gi->neighbor[v][i] > v) 
                    SetAdd(Vextension, (int)Gi->neighbor[v][i]);
            VArray[varraySize++]=v;
            ExtendSubGraph(Gi, WArray, VArray, Vextension, v, &varraySize, windowAdjList, windowRepInt, D, canonMSET);
        }
        SetFree(Vextension);
        GraphFree(Gi);
    }
    free(VArray);

    if(_windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MIN || _windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MAX)
        updateLeastFrequent(windowRepInt, canonMSET);
    MultisetFree(canonMSET);
}

void ProcessWindowRep(int *VArray, int windowRepInt) {
    // We should probably figure out a faster sort? This requires a function call for every comparison.
    int i, j;
    switch(_outputMode)
    {
        static SET* printed;
        case graphletFrequency:
            _graphletCount[windowRepInt] += _numWindowRep;
            break;
        case indexGraphlets:
            for(i=0; i<_windowSize; i++) printf("%i ",VArray[i]); 
            printf("\n%i %i\n", windowRepInt, _numWindowRep);
//          printf("\n%i\n", windowRepInt)
            for(i=0; i<_numWindowRep; i++)
            {
                for(j=0; j<_k; j++)
                    printf("%i ", _windowReps[i][j]); 
                printf("\n");
            }
            break;
        default: Abort("unknown or un-implemented outputMode");
    }
}


// This converts graphlet frequencies to concentrations or integers based on the sampling algorithm and command line arguments
void convertFrequencies(int numSamples)
{
	int i;
	if (_sampleMethod == SAMPLE_MCMC) 
	{
		if (_freqDisplayMode == count) 
		{
			for (i = 0; i < _numCanon; i++) 
			{
				_graphletCount[i] = _graphletConcentration[i] * numSamples;
			}
		}
	}
	else
	{
		if (_freqDisplayMode == concentration)
		{
			for (i = 0; i < _numCanon; i++) 
			{
				_graphletConcentration[i] = _graphletCount[i] / (double)numSamples;
			}
		}
	}
}



// This is the single-threaded BLANT function. YOU SHOULD PROBABLY NOT CALL THIS.
// Call RunBlantInThreads instead, it's the top-level entry point to call once the
// graph is finished being input---all the ways of reading input call RunBlantInThreads.
int RunBlantFromGraph(int k, int numSamples, GRAPH *G)
{
    int i,j, windowRepInt, D;
    char perm[maxK+1];
    assert(k <= G->n);
    _seed = time(0)+getpid();
    RandomSeed(_seed);
    SET *V = SetAlloc(G->n);
    TINY_GRAPH *g = TinyGraphAlloc(k);
    int varraySize = _windowSize > 0 ? _windowSize : maxK + 1;
    unsigned Varray[varraySize]; 
    InitializeConnectedComponents(G);
    //if (_sampleMethod == SAMPLE_MCMC)
	_window ? initializeMCMC(G, _windowSize, numSamples) : initializeMCMC(G, k, numSamples);
    for(i=0; i<numSamples || (_sampleFile && !_sampleFileEOF); i++)
    {
        if(_window) 
        {
            SampleGraphlet(G, V, Varray, _windowSize);
            _numWindowRep = 0; 
            if (_windowSampleMethod == WINDOW_SAMPLE_MIN || _windowSampleMethod == WINDOW_SAMPLE_MIN_D || _windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MIN)
                windowRepInt = getMaximumIntNumber(_k);
            if (_windowSampleMethod == WINDOW_SAMPLE_MAX || _windowSampleMethod == WINDOW_SAMPLE_MAX_D || _windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MAX)
                windowRepInt = -1;
            D = _k * (_k - 1) / 2;
            FindWindowRepInWindow(G, V, &windowRepInt, &D);
            if(_numWindowRep > 0)
                ProcessWindowRep(Varray, windowRepInt);
        }
        else 
        {
            SampleGraphlet(G, V, Varray, k);
            ProcessGraphlet(G, V, Varray, perm, g, k);
        }
    }
    if(_window)
    {
        for(i=0; i<_MAXnumWindowRep; i++)
            free(_windowReps[i]);
        free(_windowReps);
    }
    if (_sampleMethod == SAMPLE_MCMC && !_window)
	finalizeMCMC();
    if (_outputMode == graphletFrequency)
	convertFrequencies(numSamples);
    switch(_outputMode)
    {
	int canon;
    case indexGraphlets: case indexOrbits:
	break; // already output on-the-fly above
    case graphletFrequency:
	for(canon=0; canon<_numCanon; canon++) {
	if (_freqDisplayMode == concentration) {
	    if (SetIn(_connectedCanonicals, canon)) {
		printf("%lf ", _graphletConcentration[canon]);
		PrintCanonical(canon);
		printf("\n");
	    }
	}
	else {
	    if (SetIn(_connectedCanonicals, canon)) {
		printf("%lu ", _graphletCount[canon]);
		PrintCanonical(canon);
		printf("\n");
	    }
	}
	}
	break;
    case predictPairs:
	fprintf(stderr, "count\tnode1\tnode2\n");
	for(i=0; i < G->n-1; i++) for(j=i+1;j < G->n; j++)
	{
	    float total = (_PredictCount[i][j]+_PredictCount[j][i]);
	    if(total == 0) continue;
	    printf("%g",total);
	    if(_supportNodeNames) printf("\t%s %s\n",_nodeNames[i],_nodeNames[j]);
	    else printf("\t%d %d\n",i,j);
	}
	break;
    case outputGDV:
	for(i=0; i < G->n; i++)
	{
	    if(_supportNodeNames) printf("%s",_nodeNames[i]);
	    else printf("%d",i);
	    for(canon=0; canon < _numCanon; canon++)
		printf(" %lu", GDV(i,canon));
	    puts("");
	}
	break;
     case outputODV:
        for(i=0; i<G->n; i++) {
	    if(_supportNodeNames) printf("%s",_nodeNames[i]);
	    else printf("%d",i);
	    for(j=0; j<_numOrbits; j++) if (SetIn(_connectedCanonicals, _orbitCanonMapping[j])) {
		if (!_MCMC_UNIFORM) printf(" %lu", ODV(i,j));
		else printf(" %.12f", _doubleOrbitDegreeVector[j][i]);
		}
	    printf("\n");
	}
        break;
    default: Abort("unknown or un-implemented outputMode");
	break;
	}

#if PARANOID_ASSERTS // no point in freeing this stuff since we're about to exit; it can take significant time for large graphs.
    if(_outputMode == outputGDV) for(i=0;i<MAX_CANONICALS;i++)
	Free(_graphletDegreeVector[i]);
    if(_outputMode == outputODV) for(i=0;i<MAX_ORBITS;i++) Free(_orbitDegreeVector[i]);
	if(_outputMode == outputODV && _MCMC_UNIFORM) for(i=0;i<MAX_ORBITS;i++) Free(_doubleOrbitDegreeVector[i]);
    if(_outputMode == predictPairs) {
    }
    TinyGraphFree(g);
    SetFree(V);
    GraphFree(_G);
#endif
    if (_sampleMethod == SAMPLE_ACCEPT_REJECT)
    	fprintf(stderr,"Average number of tries per sample is %g\n", _acceptRejectTotalTries/(double)numSamples);
    return 0;
}

/*
** Fork a BLANT process and return a FILE pointer where it'll be sending stuff.
** Caller is responsible for reading all the stuff from the returned FILE pointer,
** detecting EOF on it, and fclose'ing it.
*/
FILE *ForkBlant(int k, int numSamples, GRAPH *G)
{
    int fds[2];
    pipe(fds);
    int pid = fork();
    if(pid > 0) // we are the parent
    {
	close(fds[1]); // we will not be writing to the pipe, so close it.
	return fdopen(fds[0],"r");
    }
    else if(pid == 0) // we are the child
    {
	close(fds[0]); // we will not be reading from the pipe, so close it.
	close(1); // close our usual stdout
	dup(fds[1]); // copy the write end of the pipe to fd 1.
	close(fds[1]); // close the original write end of the pipe since it's been moved to fd 1.
	RunBlantFromGraph(k, numSamples, G);
	exit(0);
	_exit(0);
	Abort("Both exit() and _exit failed???");
    }
    else
	Abort("fork failed");
}


// This is the primary entry point into BLANT, even if THREADS=1.  We assume you've already
// read the graph into G, and will do whatever is necessary to run blant with the number of
// threads specified.  Also does some sanity checking.
int RunBlantInThreads(int k, int numSamples, GRAPH *G)
{
    int c, i, j;
    assert(k == _k);
    assert(G->n >= k); // should really ensure at least one connected component has >=k nodes. TODO
    if(_outputMode == outputGDV) for(c=0;c<MAX_CANONICALS;c++)
	_graphletDegreeVector[c] = Calloc(G->n, sizeof(**_graphletDegreeVector));
    if(_outputMode == outputODV) for(i=0;i<MAX_ORBITS;i++){
	_orbitDegreeVector[i] = Calloc(G->n, sizeof(**_orbitDegreeVector));
	for(j=0;j<G->n;j++) _orbitDegreeVector[i][j]=0;
    }
    if (_outputMode == outputODV && _MCMC_UNIFORM) for(i=0;i<MAX_ORBITS;i++){
	_doubleOrbitDegreeVector[i] = Calloc(G->n, sizeof(**_doubleOrbitDegreeVector));
	for(j=0;j<G->n;j++) _doubleOrbitDegreeVector[i][j]=0.0;
    }
    if(_outputMode == predictPairs) {
	// The user uses the orbitID *relative* to the first "global" orbit listed in the orbit map,
	// so now convert that relative orbit ID to an absolute one.
	_predictOrbit1 += _orbitList[_predictOrd][0];
	_predictOrbit2 += _orbitList[_predictOrd][0];
	fprintf(stderr, "Sanity check: k %d, ordinal %d (integer %d), global orbitPair (%d,%d)\n",
	    _k, _predictOrd, _canonList[_predictOrd], _predictOrbit1, _predictOrbit2);
	TINY_GRAPH *T = TinyGraphAlloc(_k);
	for(c=0; c<_numCanon; c++) {
	    if(!SetIn(_connectedCanonicals, c)) continue;
	    int GintCanonInt = _canonList[c];
	    if(_K[GintCanonInt] != c) Fatal("_K[%d]=%d when c=%d", GintCanonInt, _K[GintCanonInt], c);
	    TinyGraphEdgesAllDelete(T);
	    Int2TinyGraph(T, GintCanonInt);
	    assert(NumReachableNodes(T, 0)==_k && T->n == _k);
	    PreProcessPredict(T, c);
#if VERBOSE
	    for(i=0;i<_k;i++) {
		for(j=0;j<_k;j++) if(_PredictMotifCount[c][i][j])
		    printf("c=%d, i=%d, j=%d, PreProc=%d\n", c,i,j, _PredictMotifCount[c][i][j]);
		puts("");
	    }
#endif
	}
    }
	
    if(_THREADS == 1)
	return RunBlantFromGraph(k, numSamples, G);

    // At this point, _THREADS must be greater than 1.
    int samplesPerThread = numSamples/_THREADS;  // will handle leftovers later if numSamples is not divisible by _THREADS

    FILE *fpThreads[_THREADS]; // these will be the pipes reading output of the parallel blants
    for(i=0;i<_THREADS;i++)
	fpThreads[i] = ForkBlant(_k, samplesPerThread, G);

    Boolean done = false;
    int lineNum = 0;
    do
    {
	#define MAX_WORDSIZE 20 // the maximum string length of a 64-bit int written in base 10
	char line[MAX_ORBITS * (MAX_WORDSIZE + 1) + 1]; // one space between words, plus a newline 
	int thread;
	for(thread=0;thread<_THREADS;thread++)	// read and then echo one line from each of the parallel instances
	{
	    fgets(line, sizeof(line), fpThreads[thread]);
	    if(feof(fpThreads[thread])) // assume if any one of them finishes, we are done.
	    {
		done = true;
		break;
	    }
	    char *nextChar = line;
	    unsigned long int count;
	    int canon, orbit, numRead;
	    switch(_outputMode)
	    {
	    case graphletFrequency:
		numRead = sscanf(line, "%lu%d", &count, &canon);
		assert(numRead == 2);
		_graphletCount[canon] += count;
		break;
	    case outputGDV:
		if(_supportNodeNames) Fatal("Oops, we don't yet support node names in multi-threaded ODV or GDV mode");
		int nodeId;
		assert(isdigit(*nextChar));
		numRead = sscanf(nextChar, "%d", &nodeId);
		assert(numRead == 1 && nodeId == lineNum);
		while(isdigit(*nextChar)) nextChar++; // read past current integer
		assert(*nextChar == ' ' || (canon == _numCanon-1 && *nextChar == '\n'));
		nextChar++;
		for(canon=0; canon < _numCanon; canon++)
		{
		    assert(isdigit(*nextChar));
		    numRead = sscanf(nextChar, "%lu", &count);
		    assert(numRead == 1);
		    GDV(lineNum,canon) += count;
		    while(isdigit(*nextChar)) nextChar++; // read past current integer
		    assert(*nextChar == ' ' || (canon == _numCanon-1 && *nextChar == '\n'));
		    nextChar++;
		}
		assert(*nextChar == '\0');
		break;
	    case outputODV:
		if(_supportNodeNames) Fatal("Oops, we don't yet support node names in multi-threaded ODV or GDV mode");
		assert(isdigit(*nextChar));
		numRead = sscanf(nextChar, "%d", &nodeId);
		assert(numRead == 1 && nodeId == lineNum);
		while(isdigit(*nextChar)) nextChar++; // read past current integer
		assert(*nextChar == ' ' || (orbit == _numOrbits-1 && *nextChar == '\n'));
		nextChar++;
		for(orbit=0; orbit < _numOrbits; orbit++)
		{
		    assert(isdigit(*nextChar));
		    numRead = sscanf(nextChar, "%lu", &count);
		    assert(numRead == 1);
		    ODV(lineNum,orbit) += count;
		    while(isdigit(*nextChar)) nextChar++; // read past current integer
		    assert(*nextChar == ' ' || (orbit == _numOrbits-1 && *nextChar == '\n'));
		    nextChar++;
		}
		assert(*nextChar == '\0');
		break;
	    case indexGraphlets: case indexOrbits:
		fputs(line, stdout);
		if(_window) 
		    while(fgets(line, sizeof(line), fpThreads[thread])) 
			fputs(line, stdout);
		break;
	    case predictPairs:
	    default:
		Abort("oops... unknown or unsupported _outputMode in RunBlantInThreads while reading child process");
		break;
	    }
	}
	lineNum++;
    } while(!done);
    for(i=0; i<_THREADS; i++)
	fclose(fpThreads[i]);  // close all the pipes

    // if numSamples is not a multiple of _THREADS, finish the leftover samples
    int leftovers = numSamples % _THREADS;
    return RunBlantFromGraph(_k, leftovers, G);
}

void BlantAddEdge(int v1, int v2)
{
    if(!_pairs) _pairs = Malloc(2*_maxEdges*sizeof(_pairs[0]));
    assert(_numEdges <= _maxEdges);
    if(_numEdges >= _maxEdges)
    {
	_maxEdges *=2;
	_pairs = Realloc(_pairs, 2*_maxEdges*sizeof(int));
    }
    _numNodes = MAX(_numNodes, v1+1); // add one since, for example, if we see a node numbered 100, numNodes is 101.
    _numNodes = MAX(_numNodes, v2+1);
    _pairs[2*_numEdges] = v1;
    _pairs[2*_numEdges+1] = v2;
    if(_pairs[2*_numEdges] == _pairs[2*_numEdges+1])
	Fatal("BlantAddEdge: edge %d (%d,%d) has equal nodes; cannot have self-loops\n", _numEdges, v1, v2);
    if(_pairs[2*_numEdges] > _pairs[2*_numEdges+1])
    {
	int tmp = _pairs[2*_numEdges];
	_pairs[2*_numEdges] = _pairs[2*_numEdges+1];
	_pairs[2*_numEdges+1] = tmp;
    }
    assert(_pairs[2*_numEdges] < _pairs[2*_numEdges+1]);
    _numEdges++;
}

int RunBlantEdgesFinished(int k, int numSamples, int numNodes, char **nodeNames)
{
    GRAPH *G = GraphFromEdgeList(_numNodes, _numEdges, _pairs, SPARSE, NULL);
    Free(_pairs);
    _nodeNames = nodeNames;
    return RunBlantInThreads(k, numSamples, G);
}

// Initialize the graph G from an edgelist; the user must allocate the pairs array
// to have 2*numEdges elements (all integers), and each entry must be between 0 and
// numNodes-1. The pairs array MUST be allocated using malloc or calloc, because
// we are going to free it right after creating G (ie., before returning to the caller.)
int RunBlantFromEdgeList(int k, int numSamples, int numNodes, int numEdges, int *pairs)
{
    assert(numNodes >= k);
    GRAPH *G = GraphFromEdgeList(numNodes, numEdges, pairs, SPARSE, NULL);
    Free(pairs);
    return RunBlantInThreads(k, numSamples, G);
}

const char const * const USAGE = \

    "USAGE: blant [-r seed] [-t threads (default=1)] [-m{outputMode}] [-d{displayMode}] {-n nSamples | -c confidence -w width} {-k k} {-w windowSize} {-s samplingMethod} {-p windowRepSamplingMethod} {graphInputFile}\n" \
    "Graph must be in one of the following formats with its extension name .\n" \
          "GML (.gml) GraphML (.xml) LGF(.lgf) CSV(.csv) LEDA(.leda) Edgelist (.el) .\n" \
    "outputMode is one of: o (ODV, the default); i (indexGraphlets); g (GDV); f (graphletFrequency).\n" \
	"-m{f}{frequencyDisplayMode} is allowed with frequencyDisplayMode being one of: i(integer or count) and d(decimal or concentration)\n" \
    "samplingMethod is one of: NBE (Node Based Expansion); EBE (Edge Based Expansion); MCMC (Markov chain Monte Carlo); RES (Lu Bressan's reservoir); AR (Accept-Reject).\n" \
    "windowRepSamplingMethod is one of: MIN (Minimizer); MAX (Maximizer); DMIN (Minimizer With Distance); DMAX (Maximizer with Distance); LFMIN (Least Frequent Minimizer); LFMAX (Least Frequent Maximizer).\n" \
	"displayMode controls how the graphlet is displayed: options are:\n"\
	"\ti (integer ordinal), d(decimal), b (binary), j (JESSE), o (ORCA).\n" \
    "At the moment, nodes must be integers numbered 0 through n-1, inclusive.\n" \
    "Duplicates and self-loops should be removed before calling BLANT.\n" \
    "k is the number of nodes in graphlets to be sampled." \
    "";

// The main program, which handles multiple threads if requested.  We simply fire off a bunch of parallel
// blant *processes* (not threads, but full processes), and simply merge all their outputs together here
// in the parent.
int main(int argc, char *argv[])
{
    int i, opt, numSamples=0;
    double confWidth = 0, confidence=0; // for confidence interval, if it's chosen

    if(argc == 1)
    {
	fprintf(stderr, "%s\n", USAGE);
	exit(1);
    }

    _THREADS = 1; 
    _k = 0;

    while((opt = getopt(argc, argv, "m:d:t:s:c:k:w:p:r:n:u")) != -1)
    {
	switch(opt)
	{
	case 'm':
	    if(_outputMode != undef) Fatal("tried to define output mode twice");
	    switch(*optarg)
	    {
	    case 'i': _outputMode = indexGraphlets; break;
	    case 'j': _outputMode = indexOrbits; break;
	    case 'p': _outputMode = predictPairs;
		char *s = optarg+1;
		_predictOrd=atoi(s);
		until(*s++==',') ;
		_predictOrbit1 = atoi(s);
		until(*s++==',') ;
		_predictOrbit2 = atoi(s);
		if(_predictOrbit2 < _predictOrbit1) {
		    int tmp = _predictOrbit2;
		    _predictOrbit2 = _predictOrbit1;
		    _predictOrbit1 = tmp;
		}
		assert(_predictOrbit1 <= _predictOrbit2);
		break;
	    case 'f': _outputMode = graphletFrequency;
		switch (*(optarg + 1))
		{
		    case 'i': _freqDisplayMode = count; break;
		    case 'd': _freqDisplayMode = concentration; break;
		    case '\0': _freqDisplayMode = freq_display_mode_undef; break;
		    default: Fatal("-mf%c: unknown frequency display mode;\n"
		    "\tmodes are i=integer(count), d=decimal(concentration)", *(optarg + 1));
		    break;
		}
	    break;
	    case 'g': _outputMode = outputGDV; break;
	    case 'o': _outputMode = outputODV; break;
	    default: Fatal("-m%c: unknown output mode;\n"
	       "\tmodes are i=indexGraphlets, j=indexOrbits, kINT,i,j=predictPairs of canonical INT columns i and j, f=graphletFrequency, g=GDV, o=ODV", *optarg);
	    break;
	    }
	    break;
	case 'd':
	    if (_displayMode != undefined) Fatal("tried to define canonical display mode twice");
	    switch(*optarg)
	    {
	    case 'b': _displayMode = binary; break;
	    case 'd': _displayMode = decimal; break;
	    case 'i': _displayMode = ordinal; break;
	    case 'j': _displayMode = jesse; break;
	    case 'o': _displayMode = orca; break;
	    default: Fatal("-d%c: unknown canonical display mode:n"
		    "\tmodes are i=integer ordinal, d=decimal, b=binary, o=orca, j=jesse", *optarg);
	    break;
	    }
	    break;
	case 't': _THREADS = atoi(optarg); assert(_THREADS>0); break;
	case 'r': _seed = atoi(optarg);
	    break;
	case 's':
	    if (_sampleMethod != -1) Fatal("Tried to define sampling method twice");
	    else if (strncmp(optarg, "NBE", 3) == 0)
		    _sampleMethod = SAMPLE_NODE_EXPANSION;
	    else if (strncmp(optarg, "FAYE", 4) == 0) 
		    _sampleMethod = SAMPLE_FAYE;
	    else if (strncmp(optarg, "EBE", 3) == 0)
		    _sampleMethod = SAMPLE_EDGE_EXPANSION;
	    else if (strncmp(optarg, "MCMC", 3) == 0)
		    _sampleMethod = SAMPLE_MCMC;
	    else if (strncmp(optarg, "RES", 3) == 0)
		    _sampleMethod = SAMPLE_RESERVOIR;
	    else if (strncmp(optarg, "AR", 2) == 0)
		    _sampleMethod = SAMPLE_ACCEPT_REJECT;
	    else
	    {
		    _sampleFileName = optarg;
		    if(strcmp(optarg,"STDIN") == 0) _sampleFile = stdin;
		    else _sampleFile = fopen(_sampleFileName, "r");
		    if(!_sampleFile)
			    Fatal("Unrecognized sampling method specified: '%s'. Options are: {NBE|EBE|MCMC|RES|FAYE|AR|{filename}}\n"
				    "If unrecognized, we try opening a file by the name '%s', but no such file exists",
				    _sampleFileName, _sampleFileName);
		    _sampleMethod = SAMPLE_FROM_FILE;
	    }
	    break;
	case 'c': confidence = atof(optarg);
	    Apology("confidence intervals not implemented yet");
	    break;
	case 'k': _k = atoi(optarg);
	    if(!(3 <= _k && _k <= 8)) Fatal("k must be between 3 and 8\n%s", USAGE);
	    break;
	case 'w':
	    _window = true;
	    _windowSize = atoi(optarg); 
	    if(_windowSize < _k) Fatal("windowSize must be at least size k\n");
	    _MAXnumWindowRep = CombinChooseDouble(_windowSize, _k);
	    _windowReps = Calloc(_MAXnumWindowRep, sizeof(int*));
	    int i;
	    for(i=0; i<_MAXnumWindowRep; i++)
		    _windowReps[i] = Calloc(_k+1, sizeof(int));
	    break;
	case 'p':
	    if (_windowSampleMethod != -1) Fatal("Tried to define window sampling method twice");
	    else if (strncmp(optarg, "DMIN", 4) == 0)
		    _windowSampleMethod = WINDOW_SAMPLE_MIN_D;
	    else if (strncmp(optarg, "DMAX", 4) == 0)
		    _windowSampleMethod = WINDOW_SAMPLE_MAX_D;
	    else if (strncmp(optarg, "MIN", 3) == 0)
		    _windowSampleMethod = WINDOW_SAMPLE_MIN;
	    else if (strncmp(optarg, "MAX", 3) == 0)
		    _windowSampleMethod = WINDOW_SAMPLE_MAX;
	    else if (strncmp(optarg, "LFMIN", 5) == 0)
		    _windowSampleMethod = WINDOW_SAMPLE_LEAST_FREQ_MIN;
	    else if (strncmp(optarg, "LFMAX", 5) == 0)
		    _windowSampleMethod = WINDOW_SAMPLE_LEAST_FREQ_MAX;
	    else
		    Fatal("Unrecognized window sampling method specified. Options are: {MAX|MIN|MAXD|MIND|LFMAX|LFMIN}\n");
	    break;
	case 'n': numSamples = atoi(optarg);
	    if(numSamples < 0) Fatal("numSamples must be non-negative\n%s", USAGE);
	    break;
	case 'u': _MCMC_UNIFORM = true;
	    break;
	default: Fatal("unknown option %c\n%s", opt, USAGE);
	}
    }
    if(_outputMode == undef) _outputMode = outputODV; // default to the same thing ORCA and Jesse us
	if (_freqDisplayMode == freq_display_mode_undef) // Default to integer(count)
		_freqDisplayMode = count;
    if (_windowSampleMethod == -1 && _window) 
       Fatal("Haven't specified window sampling method. Options are: -p{MAX|MIN|MAXD|MIND|LFMAX|LFMIN}\n");   

    if(_sampleMethod == -1) _sampleMethod = SAMPLE_MCMC;

    if(!argv[optind]) Fatal("no input graph file specified\n%s", USAGE);
    char *graphFileName = argv[optind];
    int piped = 0;
    FILE *fpGraph = readFile(argv[optind], &piped);
    if(!fpGraph) Fatal("cannot open graph input file '%s'\n", argv[optind]);
    optind++;
    assert(optind == argc);

    SetBlantDir(); // Needs to be done before reading any files in BLANT directory
    SetGlobalCanonMaps(); // needs _k to be set
    LoadMagicTable(); // needs _k to be set

    // Read it in using native Graph routine.
    _G = GraphReadEdgeList(fpGraph, SPARSE, _supportNodeNames, false);
    if(_supportNodeNames)
    {
	assert(_G->name);
	_nodeNames = _G->name;
    }
    closeFile(fpGraph, &piped);

    if(numSamples==0) numSamples = _G->numEdges * sqrt(_G->n);
    if(numSamples!=0 && confidence>0)
	Fatal("cannot specify both -s (sample size) and confidence interval");

    return RunBlantInThreads(_k, numSamples, _G);
}
