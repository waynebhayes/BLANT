#include <sys/file.h>
#include <sys/mman.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "rand48.h"
#include "blant.h"

#define PARANOID_ASSERTS 1	// turn on paranoid checking --- slows down execution by a factor of 2-3
#define USAGE "USAGE: blant [-numCores] {k} {nSamples} {graphInputFile}\n" \
    "Graph must be in edge-list format (one pair of unordered nodes on each line).\n" \
    "At the moment, nodes must be integers numbered 0 through n-1, inclusive.\n" \
    "Duplicates and self-loops should be removed before calling BLANT."

// These are mutually exclusive
#define SAMPLE_UNBIASED 0	// makes things REALLY REALLY slow.  Like 10-100 samples per second rather than a million.
#define SAMPLE_UNIF_OUTSET 1	// sample using uniform outset; about 100,000 samples per second
#define SAMPLE_CUMULATIVE 0	// Fastest, up to a million samples per second
#define MAX_TRIES 100		// max # of tries in cumulative sampling before giving up
#if (SAMPLE_UNBIASED + SAMPLE_UNIF_OUTSET + SAMPLE_CUMULATIVE) != 1
#error "must choose exactly one of the SAMPLE_XXX choices"
#endif

#define maxBk (1 << (maxK*(maxK-1)/2)) // maximum number of entries in the canon_map

// The following is the most compact way to store the permutation between a non-canonical and its canonical representative,
// when k=8: there are 8 entries, and each entry is a integer from 0 to 7, which requires 3 bits. 8*3=24 bits total.
// For simplicity we use the same 3 bits per entry, and assume 8 entries, even for k<8.  It wastes memory for k<4, but
// makes the coding much simpler.
typedef unsigned char kperm[3]; // The 24 bits are stored in 3 unsigned chars.

static unsigned int Bk, _k; // _k is the global variable storing k; Bk=actual number of entries in the canon_map for given k.

// Here's where we're lazy on saving memory, and we could do better.  We're going to allocate a static array
// that is big enough for the 256 million permutations from non-canonicals to canonicals for k=8, even if k<8.
// So we're allocating 256MBx3=768MB even if we need much less.  I figure anything less than 1GB isn't a big deal
// these days. It needs to be aligned to a page boundary since we're going to mmap the binary file into this array.
static kperm Permutations[maxBk] __attribute__ ((aligned (8192)));
// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
static short int K[maxBk] __attribute__ ((aligned (8192)));


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

// Given a TINY_GRAPH and k, return the integer ID cretaed from one triangle (upper or lower) of the adjacency matrix.
static int TinyGraph2Int(TINY_GRAPH *G, int k)
{
    int i, j, bitPos=0, Gint = 0, bit;
    
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
	    if(TinyGraphAreConnected(G,i,j))
	    {
		bit = (1 << bitPos);
		Gint |= bit;
	    }
            bitPos++;
        }
    }
    return Gint;
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
static SET *SampleGraphletUniformOutSet(SET *V, int *Varray, GRAPH *G, int k)
{
    static SET *outSet;
    if(!outSet)
       outSet = SetAlloc(G->n);  // we won't bother to free this since it's static.
    else if(G->n > outSet->n)
	SetResize(outSet, G->n);
    else
	SetEmpty(outSet);
    int edge = G->numEdges * drand48(), i, v1, v2;
    assert(V && V->n >= G->n);
    SetEmpty(V);
    int nOut = 0, outbound[G->n]; // vertices one step outside the boundary of V
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
	if(nOut == 0) // the graphlet has saturated it's connected component, return a disconnected graphette.
	{
	    while(SetIn(V, (j = G->n*drand48())))
		; // must terminate since k <= G->n
	    outbound[nOut++] = j;
	    j = 0;
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
static SET *SampleGraphletCumulative(SET *V, int *Varray, GRAPH *G, int k)
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
	while(SetIn(internal, (whichNeigh = outDegree * drand48())))
	    ; // which edge to choose among all edges leaving all nodes in V so far?
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
	    else
	    {
		// We are probably in a connected component with fewer than k nodes.
		// Test that hypothesis.
		int nodeArray[G->n], distArray[G->n];
		int sizeOfCC = GraphBFS(G, v1, G->n, nodeArray, distArray);
		assert(sizeOfCC < k);

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
#if PARANOID_ASSERTS
	assert(SetCardinality(V)==k);
#endif
	TinyGraphEdgesAllDelete(g);
	TinyGraphInducedFromGraph(g, G, Varray);
    } while(TinyGraphBFS(g, 0, k, graphetteArray, distArray) < k);

    return V;
}


// Try to mmap, and if it fails, just slurp in the file (sigh, Windoze)
void *Mmap(void *p, size_t n, int fd)
{
#if MMAP
    void *newPointer = mmap(p, n, PROT_READ, MAP_PRIVATE|MAP_FIXED, fd, 0);
    if(newPointer == MAP_FAILED)
#endif
    {
	Warning("mmap failed");
	if(read(fd, p, n) != n)
	    Fatal("cannot mmap, or cannot read the file, or both");
    }
    return p;
}

static int IntCmp(const void *a, const void *b)
{
    int *i = (int*)a, *j = (int*)b;
    return (*i)-(*j);
}

// This is the single-core version of blant. It used to be the main() function, which is why it takes (argc,argv[]).
// Now it simply gets called once for each core we use when run with multiple cores.
int blant(int argc, char *argv[])
{
    int k=atoi(argv[1]), i, j;
    _k = k;
    int numSamples = atoi(argv[2]);
    FILE *fp = fopen(argv[3], "r");
    GRAPH *G = GraphReadEdgeList(fp, true); // sparse = true
    assert(G->n >= k);
    fclose(fp);
    srand48(time(0)+getpid());

    Bk = (1 <<(k*(k-1)/2));
    char BUF[BUFSIZ];
    sprintf(BUF, CANON_DIR "/canon_list%d.txt", k);
    FILE *fp_ord=fopen(BUF, "r");
    assert(fp_ord);
    int numCanon;
    fscanf(fp_ord, "%d",&numCanon);
    int canon_list[numCanon];
    for(i=0; i<numCanon; i++) fscanf(fp_ord, "%d", &canon_list[i]);
    fclose(fp_ord);
    char perm[maxK+1];
    sprintf(BUF, CANON_DIR "/canon_map%d.bin", k);
    int Kfd = open(BUF, 0*O_RDONLY);
    sprintf(BUF, CANON_DIR "/perm_map%d.bin", k);
    int pfd = open(BUF, 0*O_RDONLY);
    short int *Kf = Mmap(K, Bk*sizeof(K[0]), Kfd);
    kperm *Pf = Mmap(Permutations, Bk*sizeof(Permutations[0]), pfd);
    assert(Kf == K);
    assert(Pf == Permutations);

    SET *V = SetAlloc(G->n);
    TINY_GRAPH *g = TinyGraphAlloc(k);
    unsigned Varray[maxK+1];
    for(i=0; i<numSamples; i++)
    {
#if SAMPLE_UNBIASED
	SampleGraphletUnbiased(V, Varray, G, k);	// REALLY REALLY SLOW
#elif SAMPLE_UNIF_OUTSET
	SampleGraphletUniformOutSet(V, Varray, G, k);
#else
#assert SAMPLE_CUMULATIVE(1)
	SampleGraphletCumulative(V, Varray, G, k); // This one is faster but less well tested and less well understood.
#endif
	// We should probably figure out a faster sort? This requires a function call for every comparison.
	qsort((void*)Varray, k, sizeof(Varray[0]), IntCmp);
	TinyGraphInducedFromGraph(g, G, Varray);
	int Gint = TinyGraph2Int(g,k);
	for(j=0;j<k;j++) perm[j]=0;
	ExtractPerm(perm, Gint);
	//printf("K[%d]=%d [%d];", Gint, K[Gint], canon_list[K[Gint]]);
	printf("%d", K[Gint]); // Note this is the ordinal of the canonical, not its bit representation
	for(j=0;j<k;j++) printf(" %d", Varray[(unsigned)perm[j]]);
	puts("");
    }
#if PARANOID_ASSERTS // no point in freeing this stuff since we're about to exit; it can take significant time for large graphs.
    TinyGraphFree(g);
    SetFree(V);
    GraphFree(G);
#endif
    return 0;
}

// The main program, which handles multiple cores if requested.  We simply fire off a bunch of parallel
// blant *processes* (not threads, but full processes), and simply merge all their outputs together here
// in the parent.
int main(int argc, char *argv[])
{
    int i, CORES;
    if(argc != 4 && argc != 5) Apology(USAGE);
    if(argc == 5) // only way to get 5 args is if the user specified "-cores" as the first argument.
    {
	CORES = atoi(argv[1]);
	if(CORES >= 0)
	    Apology("You specified 5 arguments but the first argument must be specified as -nCores\n" USAGE);
	CORES = -CORES;
	for(i=1;i<argc;i++) // nuke the numCores argument
	    argv[i]=argv[i+1];
	--argc;
    }
    else // CORES = 1;
	return blant(argc, argv);

    char cmd[BUFSIZ]; // buffer to hold the command line of the core-wise BLANT processes.
    // create the command that calls ourself CORES times, generating numSamples/CORES samples each.
    int numSamples = atoi(argv[2]);  // will handle leftovers later
    int samplesPerCore = numSamples/CORES;  // will handle leftovers later if numSamples is not divisible by CORES
    sprintf(cmd, "%s %s %d %s", argv[0], argv[1], samplesPerCore, argv[3]);

    FILE *fp[CORES]; // these will be the pipes reading output of the parallel blants
    for(i=0;i<CORES;i++)
    {
	fp[i] = popen(cmd, "r"); // fire off a single-core blant with its output to this FILE poniter.
	assert(fp[i]);
    }

    Boolean done = false;
    do
    {
	char line[BUFSIZ];
	for(i=0;i<CORES;i++)	// read and then echo one line from each of the parallel instances
	{
	    fgets(line, sizeof(line), fp[i]);
	    if(feof(fp[i])) // assume if any one of them finishes, we are done.
	    {
		done = true;
		break;
	    }
	    fputs(line, stdout);
	}
    } while(!done);
    for(i=0; i<CORES; i++)
	pclose(fp[i]);  // close all the pipes

    // if numSamples is not a multiple of CORES, finish the leftover samples
    int leftovers = numSamples % CORES;
    if(leftovers)
    {
	sprintf(cmd, "%s %s %d %s", argv[0], argv[1], leftovers, argv[3]);
	system(cmd);
    }
    return 0;
}
