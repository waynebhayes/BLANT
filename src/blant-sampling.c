#include "blant.h"
#include "blant-sampling.h"
#include "blant-utils.h"
#include "heap.h"
#include "graph.h"
#include "queue.h"
#include "multisets.h"
#include "blant-output.h"

int _sampleMethod = -1, _sampleSubmethod = -1;
FILE *_sampleFile; // if _sampleMethod is SAMPLE_FROM_FILE
char _sampleFileEOF;
Boolean _MCMC_EVERY_EDGE = false; // Should MCMC restart at each edge
int _samplesPerEdge = 0;
unsigned _MCMC_L;
unsigned long int _acceptRejectTotalTries;
GRAPH *_EDGE_COVER_G;

// Update the most recent d-graphlet to a random neighbor of it
int *MCMCGetNeighbor(int *Xcurrent, GRAPH *G)
{
    if (mcmc_d == 2)
    {
	int oldu = Xcurrent[0];
	int oldv = Xcurrent[1];
	int numTries = 0;
	while (oldu == Xcurrent[0] && oldv == Xcurrent[1]) {
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
	    if(_sampleSubmethod == SAMPLE_MCMC_EC) {
		if(!GraphAreConnected(_EDGE_COVER_G, Xcurrent[0], Xcurrent[1]) && RandomUniform() < 0.9) { // undo the attempt
		    Xcurrent[0] = oldu;
		    Xcurrent[1] = oldv;
		}
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
//   If whichCC < 0, then it's really a starting edge, where -1 means edgeList[0], -2 means edgeList[1], etc.

double SampleGraphletNodeBasedExpansion(GRAPH *G, SET *V, int *Varray, int k, int whichCC)
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
    int edge;
    if(whichCC<0){
	edge = -(whichCC+1);
	v1 = G->edgeList[2*edge];
    }
    else do {
	edge = G->numEdges * RandomUniform();
	v1 = G->edgeList[2*edge];
    } while(!SetIn(_componentSet[whichCC], v1));
    assert(edge < G->numEdges);
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
#if PARANOID_ASSERTS
	    assert(!SetIn(V, nv1)); // assertion to ensure we're in line with faye
#endif
	    SetAdd(outSet, (outbound[nOut++] = nv1));
	}
    }
    for(i=0; i < G->degree[v2]; i++)
    {
	int nv2 =  G->neighbor[v2][i];
	if(nv2 != v1 && !SetIn(outSet, nv2))
	{
#if PARANOID_ASSERTS
	    assert(!SetIn(V, nv2)); // assertion to ensure we're in line with faye
#endif
	    SetAdd(outSet, (outbound[nOut++] = nv2));
	}
    }
    for(i=2; i<k; i++)
    {
	int j;
	if(nOut == 0) // the graphlet has saturated it's connected component
	{
#if PARANOID_ASSERTS
	    assert(SetCardinality(outSet) == 0);
	    assert(SetCardinality(V) < k);
#endif
#if ALLOW_DISCONNECTED_GRAPHLETS
	    while(SetIn(V, (j = G->n*RandomUniform()))) ; // must terminate since k <= G->n
	    outbound[nOut++] = j;
	    j = 0;
#else
	    static int depth;
	    depth++;
	    // must terminate eventually as long as there's at least one connected component with >=k nodes.
	    assert(depth < MAX_TRIES); // graph is too disconnected
	    SampleGraphletNodeBasedExpansion(G, V, Varray, k, whichCC);
	    depth--;
	    // Ensure the damn thing really *is* connected.
	    TINY_GRAPH *T = TinyGraphAlloc(k);
	    TinyGraphInducedFromGraph(T, G, Varray);
#if PARANOID_ASSERTS
	    assert(NumReachableNodes(T,0) == k);
#endif
	    TinyGraphFree(T);
	    return 1.0;
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
    return 1.0;
}

// modelled after faye by Tuong Do
double SampleGraphletFaye(GRAPH *G, SET *V, int *Varray, int k, int whichCC)
{
    /* Faye: Add a visited array to keep track of nodes. Initialize to 0 */
    int visited[G->n];
    int m = 0;
    for (m = 0; m < G->n; m++) {
        visited[m] = 0;
    }
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
#if PARANOID_ASSERTS
	    assert(!SetIn(V, nv1)); // assertion to ensure we're in line with faye
#endif
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
#if PARANOID_ASSERTS
	    assert(!SetIn(V, nv2)); // assertion to ensure we're in line with faye
#endif
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
#if PARANOID_ASSERTS
	    assert(SetCardinality(outSet) == 0);
	    assert(SetCardinality(V) < k);
#endif
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
	    SampleGraphletFaye(G, V, Varray, k, whichCC);
	    depth--;
	    // Ensure the damn thing really *is* connected.
	    TINY_GRAPH *T = TinyGraphAlloc(k);
	    TinyGraphInducedFromGraph(T, G, Varray);
#if PARANOID_ASSERTS
	    assert(NumReachableNodes(T,0) == k);
#endif
	    TinyGraphFree(T);
	    return 1.0;
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
    return 1.0;
}

// Returns NULL if there are no more samples
double SampleGraphletFromFile(GRAPH *G, SET *V, int *Varray, int k)
{
    SetEmpty(V);
    int i, numRead;
    char line[BUFSIZ];
    char *s = fgets(line, sizeof(line), _sampleFile);
    if(!s){
	_sampleFileEOF = 1; // forces exit below
	return 1.0;
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
    return 1.0;
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
**   If whichCC < 0, then it's really a starting edge, where -1 means edgeList[0], -2 means edgeList[1], etc.
*/
double SampleGraphletEdgeBasedExpansion(GRAPH *G, SET *V, int *Varray, int k, int whichCC)
{
    int edge, v1, v2;
    assert(V && V->n >= G->n);
    SetEmpty(V);
    int nOut = 0;
    if(whichCC<0){
	edge = -(whichCC+1);
	v1 = G->edgeList[2*edge];
    }
    else do {
	edge = G->numEdges * RandomUniform();
	v1 = G->edgeList[2*edge];
    } while(!SetIn(_componentSet[whichCC], v1));
    assert(edge < G->numEdges);
    v2 = G->edgeList[2*edge+1];
    SetAdd(V, v1); Varray[0] = v1;
    SetAdd(V, v2); Varray[1] = v2;
    int vCount = 2;

    int outDegree = G->degree[v1] + G->degree[v2];
    static int cumulative[MAX_K];
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
	while(numTries < MAX_TRIES &&
		(whichNeigh = outDegree * RandomUniform()) >= 0 && // always true, just setting whichNeigh
		SetIn(internal, whichNeigh))
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
	    SampleGraphletEdgeBasedExpansion(G, V, Varray, k, whichCC);
	    depth--;
	    return 1.0;
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
		SampleGraphletEdgeBasedExpansion(G, V, Varray, k, whichCC);
		depth--;
		return 1.0;
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
    return 1.0;
}

// From the paper: ``Sampling Connected Induced Subgraphs Uniformly at Random''
// Xuesong Lu and Stephane Bressan, School of Computing, National University of Singapore
// {xuesong,steph}@nus.edu.sg
// 24th International Conference on Scientific and Statistical Database Management, 2012
// Note that they suggest *edge* based expansion to select the first k nodes, and then
// use reservoir sampling for the rest. But we know edge-based expansion sucks, so we'll
// start with a better starting guess, which is node-based expansion.
double SampleGraphletLuBressanReservoir(GRAPH *G, SET *V, int *Varray, int k, int whichCC)
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
	    SampleGraphletLuBressanReservoir(G, V, Varray, k, whichCC);
	    depth--;
	    return 1.0;
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
		static int graphetteArray[MAX_K], distArray[MAX_K];
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
#if PARANOID_ASSERTS
		assert(NumReachableNodes(T,0) == TinyGraphBFS(T, 0, k, graphetteArray, distArray));
#endif
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
    SampleGraphletEdgeBasedExpansion(G, V, Varray, k, whichCC);
#endif
    return 1.0;
}

/* SampleGraphletMCMC first starts with an edge in Walk L steps.
   It then walks along L = k-1 edges with MCMCGetNeighbor until it fills XLQ and XLS with L edges(their vertices are stored).
   XLQ and XLS always hold 2L vertices. They represent a window of the last L edges walked. If that window contains k
   distinct vertices, a graphlet is returned.
   This random walk predictably overcounts graphlets based on the alpha value and multipler.
   The alpha value is precomputed per graphlet type and the multiplier is based on the degree of the graphlets inside of it.
   After the algorithm is run the frequencies are normalized into concentrations.
*/

// foundGraphletCount is the expected count of the found graphlet (multiplier/_alphaList[GintOrdinal]), which needs to be returned (but must be a parameter since there's already a return value on the function)
double SampleGraphletMCMC(GRAPH *G, SET *V, int *Varray, int k, int whichCC) {
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
	if (!setup && !_MCMC_EVERY_EDGE) {
		setup = true;
		WalkLSteps(XLS, XLQ, Xcurrent, G, k, whichCC, -1);
	}
	else if (_MCMC_EVERY_EDGE && (!setup || currSamples >= _samplesPerEdge))
	{
		if(_sampleSubmethod == SAMPLE_MCMC_EC) printf("MCMC reset\n");
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
	int node, numNodes = 0, i, j;
	double multiplier = 1;
	SetEmpty(V);

	for (i = 0; i < _MCMC_L; i++) {
	    int graphletDegree = -2; // The edge between the vertices in the graphlet isn't included and is double counted
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
	    assert(multiplier > 0.0);
	}
	TinyGraphInducedFromGraph(g, G, Varray);
	Gint_type Gint = TinyGraph2Int(g, k);
	int GintOrdinal = L_K(Gint);

	assert(numNodes == k); // Ensure we are returning k nodes

	if(_sampleSubmethod == SAMPLE_MCMC_EC) {
	    int i,j;
	    for(i=0;i<k;i++) for(j=i+1;j<k;j++) {
		int u=Varray[i], v=Varray[j];
		if(GraphAreConnected(G,u,v)) GraphDisconnect(_EDGE_COVER_G,u,v);
	    }
	}
	double count = 1.0;
	if (_MCMC_L == 2) { // If _MCMC_L == 2, k = 3 and we can use the simplified overcounting formula.
	    // The over counting ratio is the alpha value only.
	    count = 1.0/(_alphaList[GintOrdinal]);
	} else {
	    // The over counting ratio is the alpha value divided by the multiplier
	    count = (double)multiplier/((double)_alphaList[GintOrdinal]);
	}
	if (_outputMode == outputODV) {
	    char perm[k];
	    memset(perm, 0, k);
	    ExtractPerm(perm, Gint);
	    for (j = 0; j < k; j++) {
		_doubleOrbitDegreeVector[_orbitList[GintOrdinal][j]][Varray[(int)perm[j]]] += count;
	    }
	} else {
	    if(count < 0) {
		Warning("count (%g) is less than 0\n", count);
	    }
	    _graphletConcentration[GintOrdinal] += count;
	}

	return count; // return the expected overcount
}

double SampleGraphletLuBressan_MCMC_MHS_without_Ooze(GRAPH *G, SET *V, int *Varray, int k) { return 1.0; } // slower
double SampleGraphletLuBressan_MCMC_MHS_with_Ooze(GRAPH *G, SET *V, int *Varray, int k) { return 1.0; } // faster!

/*
* Very slow: sample k nodes uniformly at random and throw away ones that are disconnected.
*/
double SampleGraphletAcceptReject(GRAPH *G, SET *V, int *Varray, int k)
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
    return 1.0;
}

// return how many nodes found. If you call it with startingNode == 0 then we automatically clear the visited array
static TSET _visited;
static int NumReachableNodes(TINY_GRAPH *g, int startingNode)
{
    if(startingNode == 0) TSetEmpty(_visited);
    TSetAdd(_visited,startingNode);
    unsigned int j, Varray[MAX_K], numVisited = 0;
    int numNeighbors = TSetToArray(Varray, g->A[startingNode]);
    assert(numNeighbors == g->degree[startingNode]);
    for(j=0; j<numNeighbors; j++)if(!TSetIn(_visited,Varray[j])) numVisited += NumReachableNodes(g,Varray[j]);
    return 1+numVisited;
}

// Fit SampleGraphletMCMC for windowRep implementation (alphalist and overcounting is not used here)
double SampleWindowMCMC(GRAPH *G, SET *V, int *Varray, int W, int whichCC)
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
	return 1.0;
}

/**
 * This function builds graphlets to be included into the index for the -s INDEX sampling mode. For each valid sample it takes,
 * it calls the ProcessGraphlet function to print the graphlet directly. It is called inside RunBlantFromGraph function to take
 * the given amount of samples (the amount is given with -n option) for each node in the graph as the start node
 *
 * @param G the graph
 * @param prev_nodes  the temporary set of nodes in the graphlet to build
 * @param heur_arr  the array containing the heuristic values for all nodes which is used to determine which nodes in next_step to expand to
 */
void SampleGraphletIndexAndPrint(GRAPH* G, int* prev_nodes_array, int prev_nodes_count, double *heur_arr) {
    // i, j, and neigh are just used in for loops in this function
    int i, j, neigh;
    // the tiny_graph is not used in this function. it is only used as a temporary data object as part of ProcessGraphlet (see below)
    TINY_GRAPH *g = TinyGraphAlloc(_k);

    // base case for the recursion: a k-graphlet is found, print it and return
    if (prev_nodes_count == _k) {
        // ProcessGraphlet will create the k-node induced graphlet from prev_nodes_array, and then determine if said graphlet is of a low enough multiplicity (<= multiplicity)
        // ProcessGraphlet will also check that the k nodes you passed it haven't already been printed (although, this system does not work 100% perfectly)
        // ProcessGraphlet will also print the nodes as output if the graphlet passes all checks
        ProcessGraphlet(G, NULL, prev_nodes_array, _k, g);
        return; // return here since regardless of whether ProcessGraphlet has passed or not, prev_nodes_array is already of size k so we should terminate the recursion
    }

    // Populate next_step with the following algorithm
    SET *next_step; // will contain the set of neighbors of all nodes in prev_nodes_array, excluding nodes actually in prev_nodes_array
    next_step = SetAlloc(G->n);

    // for all nodes in prev_nodes_array...
    for(i=0; i<prev_nodes_count; i++) {
        // loop through all of their neighbors and...
        for(j=0; j<G->degree[prev_nodes_array[i]]; j++) {
            neigh = G->neighbor[prev_nodes_array[i]][j];
            // if the neighbor is not in prev_nodes_array add it to the set
            if(!arrayIn(prev_nodes_array, prev_nodes_count, neigh)) {
                SetAdd(next_step, neigh); // the SET takes care of deduplication
            }
        }
    }

    // create and sort next step heur arr
    int next_step_count = SetCardinality(next_step);
    int next_step_arr[next_step_count];
#if PARANOID_ASSERTS
    assert(SetToArray(next_step_arr, next_step) == next_step_count);
#endif
    SetFree(next_step); // now that we have the next_step_arr, we no longer need the SET next_step

    // populate next_step_nwhn_arr with all nodes in next_step_arr along with their heuristic values provided in heur_arr and their names
    node_whn next_step_nwhn_arr[next_step_count]; // we only need this so that we're able to sort the array
    for (i = 0; i < next_step_count; ++i) {
        int curr_node = next_step_arr[i];
        next_step_nwhn_arr[i].node = curr_node;
        next_step_nwhn_arr[i].heur = heur_arr[curr_node];
        next_step_nwhn_arr[i].name = _nodeNames[curr_node];
    }

    int (*comp_func)(const void*, const void*);

    if (_alphabeticTieBreaking) {
        comp_func = nwhn_des_alph_comp_func;
    } else {
        comp_func = nwhn_des_rev_comp_func;
    }

    qsort((void*)next_step_nwhn_arr, next_step_count, sizeof(node_whn), comp_func); // sort by heuristic first and name second

    // Loop through neighbor nodes with Top N (-lDEGN) distinct heur values
    // If there are multiple nodes with the same heur value (which might happen with degree), we need to expand to all of them because randomly picking one to expand to would break determinism
    // If -lDEGN flag is not given, then will loop through EVERY neighbor nodes in descending order of their degree.
    int num_total_distinct_values = 0;
    double old_heur = -1; // TODO, fix this so that it's not contingent upon heuristics always being >= 0
    i = 0;
    while (i < next_step_count) {
        node_whn next_step_nwhn = next_step_nwhn_arr[i];
        double curr_heur = next_step_nwhn.heur;
        if (curr_heur != old_heur) {
            ++num_total_distinct_values;
        }
        ++i;
    }
    int num_distinct_values_to_skip = (int)(num_total_distinct_values * _topThousandth) / 1000; // algo=base
    // int num_distinct_values_to_skip = _k - prev_nodes_count - 1; // algo=stairs

    int num_distinct_values = 0;
    old_heur = -1; // TODO, fix this so that it's not contingent upon heuristics not being -1
    i = 0;
    while (i < next_step_count) {
        node_whn next_step_nwhn = next_step_nwhn_arr[i];
        double curr_heur = next_step_nwhn.heur;
        if (curr_heur != old_heur) {
            ++num_distinct_values;
        }
        old_heur = curr_heur;

        // continue for the heurs we skip
        if (num_distinct_values <= num_distinct_values_to_skip) {
            ++i;
            continue;
        }

        // break once we've gotten enough distinct heur values
        if (_numWindowRepLimit != 0 && num_distinct_values - num_distinct_values_to_skip > _numWindowRepLimit) {
            break;
        }

        // perform the standard DFS step of set next, recurse with size + 1, and then unset next
        prev_nodes_array[prev_nodes_count] = next_step_nwhn.node;
        SampleGraphletIndexAndPrint(G, prev_nodes_array, prev_nodes_count + 1, heur_arr);
        // the "unset" step is commented out for efficiency since it's not actually necessary, but the comment improves readability
        // prev_nodes_array[prev_nodes_count] = 0;
        ++i;
    }
}


// if cc == G->n, then we choose it randomly. Otherwise use cc given.
double SampleGraphlet(GRAPH *G, SET *V, unsigned Varray[], int k, int cc) {
    double randomComponent = RandomUniform();
    double overcount = 1.0; // this will be returned
    if(cc == G->n) for(cc=0; cc<_numConnectedComponents;cc++)
	if(_cumulativeProb[cc] > randomComponent)
	    break;
    //printf("choosing from CC %d\n", cc);
    switch (_sampleMethod) {
    case SAMPLE_ACCEPT_REJECT:
	SampleGraphletAcceptReject(G, V, Varray, k);	// REALLY REALLY SLOW and doesn't need to use cc
	break;
    case SAMPLE_NODE_EXPANSION:
	SampleGraphletNodeBasedExpansion(G, V, Varray, k, cc);
	break;
    case SAMPLE_FAYE:
	SampleGraphletFaye(G, V, Varray, k, cc);
	break;
    case SAMPLE_RESERVOIR:
	SampleGraphletLuBressanReservoir(G, V, Varray, k, cc); // pretty slow but not as bad as unbiased
	break;
    case SAMPLE_EDGE_EXPANSION:
	SampleGraphletEdgeBasedExpansion(G, V, Varray, k, cc); // Faster than NBE but less well tested and understood.
	break;
    case SAMPLE_MCMC:
        if(!_window) {
	    overcount = SampleGraphletMCMC(G, V, Varray, k, cc);
        } else {
            SampleWindowMCMC(G, V, Varray, k, cc);
        }
	break;
    case SAMPLE_FROM_FILE:
	SampleGraphletFromFile(G, V, Varray, k);
	break;
    case -1:
	Fatal("Please specify a sampling method using the '-s' option");
	break;
    default:
	Fatal("unknown sampling method");
	break;
    }

    return overcount;
}
