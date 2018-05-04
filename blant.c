#include <sys/file.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "rand48.h"
#include "heap.h"
#include "blant.h"
#include "queue.h"

#define PARANOID_ASSERTS 1	// turn on paranoid checking --- slows down execution by a factor of 2-3

// Enable the code that uses C++ to parse input files?
#define SHAWN_AND_ZICAN 0
static int *_pairs, _numNodes, _numEdges, _maxEdges=1024;
char **_nodeNames;

// Below are the sampling methods; pick one on the last line
#define SAMPLE_UNBIASED 1	// makes things REALLY REALLY slow.  Like 10-100 samples per second rather than a million.
#define SAMPLE_NODE_EXPANSION 2	// sample using uniform node expansion; about 100,000 samples per second
#define SAMPLE_EDGE_EXPANSION 3	// Fastest, up to a million samples per second
#define SAMPLE_RESERVOIR 4	// Lu Bressan's reservoir sampler, reasonably but not entirely unbiased.
#define SAMPLE_MCMC 5 // MCMC Algorithm estimates graphlet frequency with a random walk
#ifndef RESERVOIR_MULTIPLIER
// this*k is the number of steps in the Reservoir walk. 8 seems to work best, empirically.
#define RESERVOIR_MULTIPLIER 8
#endif
#define SAMPLE_METHOD SAMPLE_MCMC

#define MAX_TRIES 100		// max # of tries in cumulative sampling before giving up

#define ALLOW_DISCONNECTED_GRAPHLETS 0

// The following is the most compact way to store the permutation between a non-canonical and its canonical representative,
// when k=8: there are 8 entries, and each entry is a integer from 0 to 7, which requires 3 bits. 8*3=24 bits total.
// For simplicity we use the same 3 bits per entry, and assume 8 entries, even for k<8.  It wastes memory for k<4, but
// makes the coding much simpler.
typedef unsigned char kperm[3]; // The 24 bits are stored in 3 unsigned chars.

static unsigned int _Bk, _k; // _k is the global variable storing k; _Bk=actual number of entries in the canon_map for given k.
static int _numCanon, _canonList[MAX_CANONICALS];
static int _numOrbits, _orbitList[MAX_CANONICALS][maxK]; // Jens: this may not be the array we need, but something like this...

enum OutputMode {undef, indexGraphlets, graphletFrequency, outputODV, outputGDV};
static enum OutputMode _outputMode = undef;
static unsigned long int _graphletCount[MAX_CANONICALS];

// A bit counter-intuitive: we need to allocate this many vectors each of length [_numNodes],
// and then the degree for node v, graphlet/orbit g is _degreeVector[g][v], NOT [v][g].
// We do this simply because we know the length of MAX_CANONICALS so we pre-know the length of
// the first dimension, otherwise we'd need to get more funky with the pointer allocation.
// Only one of these actually get allocated, depending upon outputMode.
static unsigned long int *_graphletDegreeVector[MAX_CANONICALS];
static unsigned long int    *_orbitDegreeVector[MAX_ORBITS];
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
static short int _K[maxBk] __attribute__ ((aligned (8192)));

static unsigned L; // walk length for MCMC algorithm
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


// return how many nodes found. If you call it with startingNode == 0 then we automatically clear the visited array
static TSET _visited;
static int NumReachableNodes(TINY_GRAPH *g, int startingNode)
{
    if(startingNode == 0) TSetEmpty(_visited);
    TSetAdd(_visited,startingNode);
    int j, Varray[maxK], numVisited = 0;
    int numNeighbors = TSetToArray(Varray, g->A[startingNode]);
    assert(numNeighbors == g->degree[startingNode]);
    for(j=0; j<numNeighbors; j++)if(!TSetIn(_visited,Varray[j])) numVisited += NumReachableNodes(g,Varray[j]);
    return 1+numVisited;
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

static SET *SampleGraphletNodeBasedExpansion(SET *V, int *Varray, GRAPH *G, int k)
{
    static SET *outSet;
    static int numIsolatedNodes;
    if(!outSet)
       outSet = SetAlloc(G->n);  // we won't bother to free this since it's static.
    else if(G->n > outSet->n)
	SetResize(outSet, G->n);
    else
	SetEmpty(outSet);
    int v1, v2, i;
    int nOut = 0, outbound[G->n]; // vertices one step outside the boundary of V
    assert(V && V->n >= G->n);
    SetEmpty(V);
    int edge = G->numEdges * drand48();
    v1 = G->edgeList[2*edge];
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
    for(i=2; i<k; i++)
    {
	int j;
	if(nOut == 0) // the graphlet has saturated it's connected component
	{
#if ALLOW_DISCONNECTED_GRAPHLETS
	    while(SetIn(V, (j = G->n*drand48())))
		; // must terminate since k <= G->n
	    outbound[nOut++] = j;
	    j = 0;
#else
	    static int depth;
	    depth++;
	    // must terminate eventually as long as there's at least one connected component with >=k nodes.
	    assert(depth < MAX_TRIES); // graph is too disconnected
	    V = SampleGraphletNodeBasedExpansion(V, Varray, G, k);
	    depth--;
	    return V;
#endif
	}
	else
	    j = nOut * drand48();
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
static SET *SampleGraphletEdgeBasedExpansion(SET *V, int *Varray, GRAPH *G, int k)
{
    int edge = G->numEdges * drand48(), v1, v2;
    assert(V && V->n >= G->n);
    SetEmpty(V);
    int nOut = 0;
    v1 = G->edgeList[2*edge];
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
	int i, whichNeigh;
	while(numTries < MAX_TRIES && SetIn(internal, (whichNeigh = outDegree * drand48())))
	    ++numTries; // which edge to choose among all edges leaving all nodes in V so far?
	if(numTries >= MAX_TRIES) {
#if ALLOW_DISCONNECTED_GRAPHLETS
		    // get a new node outside this connected component.
		    // Note this will return a disconnected graphlet.
		    while(SetIn(V, (newNode = G->n*drand48())))
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
		    V = SampleGraphletEdgeBasedExpansion(V, Varray, G, k);
		    depth--;
		    return V;
#endif
	}
	for(i=0; cumulative[i] <= whichNeigh; i++)
	    ; // figure out whose neighbor it is
	int localNeigh = whichNeigh-(cumulative[i]-G->degree[Varray[i]]); // which neighbor of node i?
	int newNode = G->neighbor[Varray[i]][localNeigh];
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
		while(SetIn(V, (newNode = G->n*drand48())))
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
		V = SampleGraphletEdgeBasedExpansion(V, Varray, G, k);
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
static SET *SampleGraphletLuBressanReservoir(SET *V, int *Varray, GRAPH *G, int k)
{
    // Start by getting the first k nodes using a previous method. Once you figure out which is
    // better, it's probably best to share variables so you don't have to recompute the outset here.
#if 1  // the following is copied almost verbatim from NodeEdgeExpansion, just changing for loop to while loop.
    static SET *outSet;
    static int numIsolatedNodes;
    if(!outSet)
       outSet = SetAlloc(G->n);  // we won't bother to free this since it's static.
    else if(G->n > outSet->n)
	SetResize(outSet, G->n);
    else
	SetEmpty(outSet);
    int v1, v2, i;
    int nOut = 0, outbound[G->n]; // vertices one step outside the boundary of V
    assert(V && V->n >= G->n);
    SetEmpty(V);
    int edge = G->numEdges * drand48();
    v1 = G->edgeList[2*edge];
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
		while(SetIn(V, (v1 = G->n*drand48())))
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
	    V = SampleGraphletLuBressanReservoir(V, Varray, G, k);
	    depth--;
	    return V;
#endif
	}
	else
	{
	    candidate = nOut * drand48();
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
	    double reservoir_alpha = drand48();
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
		int memberToDelete = k*drand48();
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
    SampleGraphletEdgeBasedExpansion(V, Varray, G, k);
#endif
}

// MCMC getNeighbor
int *MCMCGetNeighbor(int *Xcurrent, GRAPH *G)
{
	//*modularize
	static SET *outSet;
    static int numIsolatedNodes;
    if(!outSet)
    	outSet = SetAlloc(G->n);
    else if(G->n > outSet->n)
		SetResize(outSet, G->n);
    else
		SetEmpty(outSet);

	int nOut = 0, outbound[G->n], p, j, i;
	
	// if (mcmc_d != 2)
	// {
	// 	j =  X->n * drand48();
	// 	while(true)
	// 	{
	// 		int vx = Varray[j];
	// 		//*modularize
	// 		for(i=0; i < G->degree[v1]; i++)
	// 		{
	// 			int nv1 =  G->neighbor[v1][i];
	// 			if(nv1 != vx) SetAdd(outSet, (outbound[nOut++] = nv1));
	// 		}

	// 		j = nOut * drand48();
	// 		vj = outbound[j];
	// 		SetDelete(X, vx);
	// 		SetAdd(X, vj); // also add to array, does local set/arr need to be updated
			
	// 		//vj;
	// 		// Xnext = swap vx & vj
	// 		// if Xnext is connected, break
	// 		break;
	// 	}
	// }
	
	if (mcmc_d == 2)
	{
		// assuming d = 2 (2 node graphlet)
		// randomly create a probability
		int oldu = Xcurrent[0];
		int oldv = Xcurrent[1];
		while (oldu == Xcurrent[0] && oldv == Xcurrent[1]) {
			p = drand48();
			// if 0 < p < 1, p < deg(u) + deg(v) then
			if (p < G->degree[Xcurrent[0]]/(G->degree[Xcurrent[0]] + G->degree[Xcurrent[1]])) {
				// select randomly from Neigh(u) and swap
				int neighbor = (int) (G->degree[Xcurrent[0]] * drand48());
				Xcurrent[1] = G->neighbor[Xcurrent[0]][neighbor];
			}
			else {
				// select randomly from Neigh(v) and swap
				int neighbor = (int) (G->degree[Xcurrent[1]] * drand48());
				Xcurrent[0] = G->neighbor[Xcurrent[1]][neighbor];
			}
		}
	}
	else Fatal("Not implemented. Set d to 2");
	return Xcurrent;
}

// WalkLSteps
SET *WalkLSteps(int *Varray, SET *UnionNodes, QUEUE *XL, GRAPH *G, int k)
{
	static int X[mcmc_d];
	int numNodes = 0;
	int edge = G->numEdges * drand48();
    int v1 = G->edgeList[2*edge];
    int v2 = G->edgeList[2*edge+1];
	Varray[0] = v1;
	Varray[1] = v2;
    SetAdd(UnionNodes, v1); QueuePut(XL, (foint) v1);
    SetAdd(UnionNodes, v2); QueuePut(XL, (foint) v2);
	X[numNodes++] = v1;
	X[numNodes++] = v2;
	while (SetCardinality(UnionNodes) < k) {
		MCMCGetNeighbor(X, G);
		if (!SetIn(UnionNodes, X[0])) Varray[numNodes++] = X[0];
		if (!SetIn(UnionNodes, X[1])) Varray[numNodes++] = X[1];
		SetAdd(UnionNodes, X[0]); QueuePut(XL, (foint) X[0]);
		SetAdd(UnionNodes, X[1]); QueuePut(XL, (foint) X[1]);

	}
	return UnionNodes;
}

void FillSetFromQueue(SET* V, QUEUE *Q)
{
	int num;
	for (int i = 0; i < Q->length; i++) {
		int num = (Q->queue[(Q->front + i) % Q->maxSize]).i;
		if (!SetIn(V, num)) {
			SetAdd(V, num);
		}
	}
}

// MCMC sampleGraphletMCMC
static SET *SampleGraphletMCMC(SET *V, int *Varray, GRAPH *G, int k)
{
	static Boolean setup = false;
	static SET *UnionNodes;
	static QUEUE *XL;

	if (!XL) {
		XL = QueueAlloc(_k*2);
	}

	if (!setup) {
		return WalkLSteps(Varray, V, XL, G, k);
	}
	else {
		int v1, v2;
		do {
			QueueGet(XL); QueueGet(XL); //Pop off the oldest graphlet from the queue
			SetEmpty(V);
			FillSetFromQueue(V, XL);
		}
		while(SetCardinality(V) >= k);

		while (SetCardinality(V) < k)
		{
			v1 = (XL->queue[(XL->front + (XL->length - 1)) % XL->maxSize]).i;
			v2 = (XL->queue[(XL->front + (XL->length)) % XL->maxSize]).i;
			SetAdd(V, v1); QueuePut(XL, (foint) v1);
			SetAdd(V, v2); QueuePut(XL, (foint) v2);
		}
		SetEmpty(V);

		int num, numNodes = 0;
		for (int i = 0; i < XL->length; i++) {
			int num = (XL->queue[(XL->front + i) % XL->maxSize]).i;
			if (!SetIn(V, num)) {
				Varray[numNodes++] = num;
				SetAdd(V, num);
			}
		}
		return MCMCGETNeighbor(G, Varray);
	}

}

void initializeMCMC() {
	L = _k - mcmc_d  + 1;
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

static SET *SampleGraphletUnbiased(SET *V, int *Varray, GRAPH *G, int k)
{
    int arrayV[k], i;
    int nodeArray[G->n], distArray[G->n];
    TINY_GRAPH *g = TinyGraphAlloc(k);
    int graphetteArray[k];

    do
    {
	SetEmpty(V);
	// select k nodes uniformly at random from G without regard to connectivity
	for(i=0; i<k; i++)
	{
	    do
		Varray[i] = G->n * drand48();
	    while(SetIn(V, Varray[i]));
	    SetAdd(V, Varray[i]);
	}
	TinyGraphEdgesAllDelete(g);
	TinyGraphInducedFromGraph(g, G, Varray);
#if PARANOID_ASSERTS
	assert(SetCardinality(V)==k);
	assert(NumReachableNodes(g,0) == TinyGraphBFS(g, 0, k, graphetteArray, distArray));
#endif
    } while(NumReachableNodes(g, 0) < k);
    return V;
}


static int IntCmp(const void *a, const void *b)
{
    int *i = (int*)a, *j = (int*)b;
    return (*i)-(*j);
}


// Assuming the global variable _k is set properly, go read in and/or mmap the big global
// arrays related to canonical mappings and permutations.
void SetGlobalCanonMaps(void)
{
    assert(3 <= _k && _k <= 8);
    _Bk = (1 <<(_k*(_k-1)/2));
    char BUF[BUFSIZ];
    int i;
    _numCanon = canonListPopulate(BUF, _canonList, _k);
    // Jens: this is where you'd insert a _numOrbits = orbitListPopulate(...) function;
    // put the actual function in libblant.[ch]
    mapCanonMap(BUF, _K, _k);
    sprintf(BUF, CANON_DIR "/perm_map%d.bin", _k);
    int pfd = open(BUF, 0*O_RDONLY);
    kperm *Pf = Mmap(Permutations, _Bk*sizeof(Permutations[0]), pfd);
    assert(Pf == Permutations);
}

void SampleGraphlet(GRAPH *G, SET *V, unsigned Varray[], int k) {
#if SAMPLE_METHOD == SAMPLE_UNBIASED
	SampleGraphletUnbiased(V, Varray, G, k);	// REALLY REALLY SLOW
#elif SAMPLE_METHOD == SAMPLE_NODE_EXPANSION
	SampleGraphletNodeBasedExpansion(V, Varray, G, k);
#elif SAMPLE_METHOD == SAMPLE_RESERVOIR
	SampleGraphletLuBressanReservoir(V, Varray, G, k); // pretty slow but not as bad as unbiased
#elif SAMPLE_METHOD == SAMPLE_EDGE_EXPANSION
	SampleGraphletEdgeBasedExpansion(V, Varray, G, k); // This one is faster but less well tested and less well understood.
#elif SAMPLE_METHOD == SAMPLE_MCMC
	SampleGraphletMCMC(V, Varray, G, k);
#else
#error unknown sampling method
#endif
}

void ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], char perm[], TINY_GRAPH *g, int k) {
	// We should probably figure out a faster sort? This requires a function call for every comparison.
	qsort((void*)Varray, k, sizeof(Varray[0]), IntCmp);
	TinyGraphInducedFromGraph(g, G, Varray);
	int Gint = TinyGraph2Int(g,k), j, GintCanon=_K[Gint];
#if PARANOID_ASSERTS
	assert(0 <= GintCanon && GintCanon < _numCanon);
#endif
	switch(_outputMode)
	{
	case graphletFrequency:
	    ++_graphletCount[GintCanon];
	    break;
	case indexGraphlets:
	    memset(perm, 0, k);
	    ExtractPerm(perm, Gint);
	    printf("%d", GintCanon); // Note this is the ordinal of the canonical, not its bit representation
#if SHAWN_AND_ZICAN
	    for(j=0;j<k;j++) printf(" %s", _nodeNames[Varray[(int)perm[j]]]);
#else
	    for(j=0;j<k;j++) printf(" %d", Varray[(int)perm[j]]);
#endif
	    puts("");
	    break;
	case outputGDV:
	    for(j=0;j<k;j++) ++GDV(Varray[j], GintCanon);
	    break;
	case outputODV: // Jens, here...
	    Abort("Sorry, ODV is not implemented yet");
	    for(j=0;j<k;j++) ++ODV(Varray[j], _orbitList[GintCanon][j]);
	    break;
		
	default: Abort("unknown or un-implemented outputMode");
	    break;
	}
}

// This is the single-threaded BLANT function. YOU SHOULD PROBABLY NOT CALL THIS.
// Call RunBlantInThreads instead, it's the top-level entry point to call once the
// graph is finished being input---all the ways of reading input call RunBlantInThreads.
int RunBlantFromGraph(int k, int numSamples, GRAPH *G)
{
    int i;
    char perm[maxK+1];
    assert(k <= G->n);
    srand48(time(0)+getpid());
    SET *V = SetAlloc(G->n);
    TINY_GRAPH *g = TinyGraphAlloc(k);
    unsigned Varray[maxK+1];
#if SAMPLE_METHOD == SAMPLE_MCMC
	initializeMCMC();
#endif
    for(i=0; i<numSamples; i++)
    {
		SampleGraphlet(G, V, Varray, k);
		ProcessGraphlet(G, V, Varray, perm, g, k);
    }
    switch(_outputMode)
    {
	int canon;
    case indexGraphlets: 
#if SAMPLE_METHOD == SAMPLE_MCMC
	Warning("Sampling method MCMC overcounts graphlets by varying amounts.");
#endif
	break; // already output on-the-fly above
    case graphletFrequency:
	for(canon=0; canon<_numCanon; canon++)
	    printf("%lu %d\n", _graphletCount[canon], canon);
	break;
    case outputGDV:
	for(i=0; i < G->n; i++)
	{
	    printf("%lu", GDV(i,0));
	    for(canon=1; canon < _numCanon; canon++)
		printf(" %lu", GDV(i,canon));
	    puts("");
	}
	break;
    default: Abort("unknown or un-implemented outputMode");
	break;
    }
#if PARANOID_ASSERTS // no point in freeing this stuff since we're about to exit; it can take significant time for large graphs.
    TinyGraphFree(g);
    SetFree(V);
    GraphFree(G);
    if(_outputMode == outputGDV) for(i=0;i<MAX_CANONICALS;i++)
	Free(_graphletDegreeVector[i]);
    if(_outputMode == outputODV) for(i=0;i<MAX_ORBITS;i++)
	Free(_orbitDegreeVector[i]);
#endif
    return 0;
}

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
    }
    else
	Abort("fork failed");
}


// This is the primary entry point into BLANT, even if THREADS=1.  We assume you've already
// read the graph into G, and will do whatever is necessary to run blant with the number of
// threads specified.  Also does some sanity checking.
int RunBlantInThreads(int k, int numSamples, GRAPH *G)
{
    int i;
    assert(k == _k);
    assert(G->n >= k); // should really ensure at least one connected component has >=k nodes. TODO
    if(_outputMode == outputGDV) for(i=0;i<MAX_CANONICALS;i++)
	_graphletDegreeVector[i] = Calloc(G->n, sizeof(**_graphletDegreeVector));
    if(_outputMode == outputODV) for(i=0;i<MAX_ORBITS;i++)
	_orbitDegreeVector[i] = Calloc(G->n, sizeof(**_orbitDegreeVector));

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
	    int canon, numRead;
	    switch(_outputMode)
	    {
	    case graphletFrequency:
		numRead = sscanf(line, "%lu%d", &count, &canon);
		assert(numRead == 2);
		_graphletCount[canon] += count;
		break;
	    case outputGDV:
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
	    case indexGraphlets:
		fputs(line, stdout);
		break;
	    default:
		Abort("oops... unknown _outputMode in RunBlantInThreads while reading child process");
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
    GRAPH *G = GraphFromEdgeList(_numNodes, _numEdges, _pairs, true); // sparse = true
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
    GRAPH *G = GraphFromEdgeList(numNodes, numEdges, pairs, true); // sparse = true
    Free(pairs);
    return RunBlantInThreads(k, numSamples, G);
}


const char const * const USAGE = \
    "USAGE: blant [-t threads (default=1)] [-m{outputMode}] {-s nSamples | -c confidence -w width} {-k k} {graphInputFile}\n" \
    "Graph must be in one of the following formats with its extension name .\n" \
          "GML (.gml) GraphML (.xml) LGF(.lgf) CSV(.csv) LEDA(.leda) Edgelist (.el) .\n" \
    "outputMode is one of: o (ODV, the default); i (indexGraphlets); g (GDV); f (graphletFrequency).\n" \
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

    while((opt = getopt(argc, argv, "m:t:s:c:w:k:")) != -1)
    {
	switch(opt)
	{
	case 'm':
	    if(_outputMode != undef) Fatal("tried to define output mode twice");
	    switch(*optarg)
	    {
	    case 'i': _outputMode = indexGraphlets; break;
	    case 'f': _outputMode = graphletFrequency; break;
	    case 'g': _outputMode = outputGDV; break;
	    case 'o': _outputMode = outputODV; break;
	    default: Fatal("-m%c: unknown output mode;\n"
		   "\tmodes are i=indexGraphlets, f=graphletFrequency, g=GDV, o=ODV", *optarg);
		break;
	    }
	    break;
	case 't': _THREADS = atoi(optarg); assert(_THREADS>0); break;
	case 's': numSamples = atoi(optarg);
	    if(numSamples < 0) Fatal("numSamples must be non-negative\n%s", USAGE);
	    break;
	case 'c': confidence = atof(optarg);
	    if(confidence <= 0) Fatal("-c argument (confidence of confidence interval) must be positive\n%s", USAGE);
	    break;
	case 'w': confWidth = atof(optarg);
	    if(confWidth <= 0) Fatal("-w argument (width of confidence interval) must be positive\n%s", USAGE);
	    break;
	case 'k': _k = atoi(optarg);
	    if(!(3 <= _k && _k <= 8)) Fatal("k must be between 3 and 8\n%s", USAGE);
	    break;
	default: Fatal("unknown option %c\n%s", opt, USAGE);
	}
    }
    if(_outputMode == undef) _outputMode = outputODV; // default to the same thing ORCA and Jesse use

    if(confidence>0) assert(confWidth>0);
    if(confWidth>0) assert(confidence>0);
    if(confidence>0) Apology("confidence intervals not implemented yet");

    if(numSamples!=0 && confidence>0)
	Fatal("cannot specify both -s (sample size) and confidence interval (-w, -c) pair");

    if(!argv[optind]) Fatal("no input graph file specified\n%s", USAGE);
    char *graphFileName = argv[optind];
    FILE *fpGraph = fopen(argv[optind], "r");
    if(!fpGraph) Fatal("cannot open graph input file '%s'\n", argv[optind]);
    optind++;
    assert(optind == argc);

    SetGlobalCanonMaps(); // needs _k to be set

#if SHAWN_AND_ZICAN
#ifdef CPP_CALLS_C  // false by default
    while(!feof(fpGraph))
    {
	static int line;
	int v1, v2;
	++line;
	if(fscanf(fpGraph, "%d%d ", &v1, &v2) != 2)
	    Fatal("can't find 2 ints on line %d\n", line);
	BlantAddEdge(v1, v2);
    }
    fclose(fpGraph);
#else // Shawn + Zican see here:
    fclose(fpGraph);
    _nodeNames = convertToEL(graphFileName);
    assert(_numNodes > 0);
    assert(_nodeNames && _nodeNames[0]);
    //assert(!_nodeNames[_numNodes]);
#if 1
    for(i=0; i < _numNodes; i++)
	printf("nodeName[%d]=%s\n", i, _nodeNames[i]);
    exit(0);
#endif
    // call clean maybe?
#endif
    return RunBlantEdgesFinished(_k, numSamples, _numNodes, NULL);
#else
    // Read it in using native Graph routine.
    GRAPH *G = GraphReadEdgeList(fpGraph, true); // sparse = true
    fclose(fpGraph);
    return RunBlantInThreads(_k, numSamples, G);
#endif
}
