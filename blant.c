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

#define PARANOID 1	// turn on paranoid checking --- slows down execution by a factor of 2-3

#define maxBk (1 << (maxK*(maxK-1)/2)) // maximum number of entries in the canon_map

// The following is the most compact way to store the permutation between a non-canonical and its canonical representative,
// when k=8: there are 8 entries, and each entry is a integer from 0 to 7, which requires 3 bits. 8*3=24 bits total.
// For simplicity we use the same 3 bits per entry, and assume 8 entries, even for k<8.  It wastes a *bit* of memory but
// makes the coding much simpler.
typedef unsigned char kperm[3]; // The 24 bits are stored in 3 unsigned chars.

static unsigned int Bk, _k; // _k is the global variable storing k; Bk=actual number of entries in the canon_map for given k.

// Here's where we're lazy on saving memory, and we could do better.  We're going to allocate a static array
// that is big enough for the 256 million permutations from non-canonicals to canonicals for k=8, even if k<8.
// So we're allocating 800MB even if we need much less.  I figure anything less than 1GB isn't a big deal these days.
// It needs to be aligned to a page boundary since we're going to mmap the binary file into this array.
static kperm Permutations[maxBk] __attribute__ ((aligned (8192)));
// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
static short int K[maxBk] __attribute__ ((aligned (8192)));

// you provide a permutation array, we fill it with the permutation extracted from the compressed Permutation mapping.
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
    
#if LOWER_TRIANGLE
    for(i=k-1;i>0;i--)
    {
        for(j=i-1;j>=0;j--)
#else   // UPPER_TRIANGLE
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

// Given the big graph G and an integer k, return a k-graphlet from G.
// Caller is responsible for allocating the set V and its array Varray.
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
	if(nOut == 0) // the graphlet has saturated it's connected component
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
#if PARANOID
    assert(SetCardinality(V) == k);
#endif
    return V;
}


/* Basic idea: in order to avoid actually building the outSet as we do
** above, we instead simply keep a *estimated* count C[v] of separate
** outSets of each node v in the accumulating graphlet. Then we pick
** an integer in the range [0,..sum(C)), go find out what node that is,
** and if it turns out to be a node already in V, then we simply try again.
*/
static SET *SampleGraphletCumulativeDist(SET *V, int *Varray, GRAPH *G, int k)
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

    while(vCount < k)
    {
	int i, whichNeigh = outDegree * drand48(); // which of the list of all neighbors of nodes in V so far?
	for(i=0; cumulative[i] <= whichNeigh; i++)
	    ;
	assert(i < vCount);
	int localNeigh = whichNeigh-(cumulative[i]-G->degree[Varray[i]]);
	assert(0 <= localNeigh && localNeigh < G->degree[Varray[i]]);
	int newNode = G->neighbor[Varray[i]][localNeigh];
	assert(0 <= newNode && newNode < G->n);
	if(!SetIn(V, newNode)) // yipee! Found a new node!
	{
	    SetAdd(V, newNode);
	    cumulative[vCount] = G->degree[newNode] + cumulative[vCount-1];
	    Varray[vCount++] = newNode;
#if PARANOID
	    assert(SetCardinality(V) == vCount);
#endif
	    outDegree += G->degree[newNode];
	    assert(outDegree == cumulative[vCount-1]);
	}
    }
#if PARANOID
    assert(SetCardinality(V) == k);
#endif
    assert(vCount == k);
    return V;
}

static SET *SampleGraphletUnbiasedMaybe(GRAPH *G, int k)
{
    SET *V = SetAlloc(G->n);
    int arrayV[k], i;
    int nodeArray[G->n], distArray[G->n];
    arrayV[0] = G->n * drand48();
#define MAX_TRIES 100
    int numBFS, tries=0;
    
    while((numBFS = GraphBFS(G, arrayV[0], k-1, nodeArray, distArray)) < k) {
	if(++tries > MAX_TRIES)
	    Fatal("can't find enough nodes in BFS after %d tries", MAX_TRIES);
	arrayV[0] = G->n * drand48(); // try another seed
	// NOTE: can speed this up by initializing an array with *which* nodes actually have k-1 nodes within dist k-1,
	// and then only picking from that array.  We could also fail immediately at start time if we find no such nodes.
    }
    SetAdd(V, arrayV[0]);
    
    // Now pick the rest from the nodeArray
    for(i=1; i<k; i++)
    {
	int element = 1+(numBFS-1)*drand48();
	arrayV[i] = nodeArray[element];
	assert(numBFS>0);
	nodeArray[element] = nodeArray[--numBFS];
	SetAdd(V, arrayV[i]);
    }
#if PARANOID
    assert(SetCardinality(V) == k);
#endif
    // Now check that it's connected
    GRAPH *graphette = GraphInduced(NULL, G, V);
    if(GraphBFS(graphette, 0, k, nodeArray, distArray) < k)
    {
	GraphFree(graphette);
	SetFree(V);
	return SampleGraphletUnbiasedMaybe(G, k);
    }
    GraphFree(graphette);
    return V;
}


// Try to mmap, and if it fails, just slurp in the file (sigh, Windoze)
void *Mmap(void *p, size_t n, int fd)
{
#if MMAP
    void *newPointer = mmap(p, n, PROT_READ, MAP_PRIVATE|MAP_FIXED, fd, 0);
    if(newPointer == MAP_FAILED)
#endif
	if(read(fd, p, n) != n)
	    Fatal("cannot mmap, or cannot read the file, or both");
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
#if PARANOID
	SampleGraphletUniformOutSet(V, Varray, G, k);
#else
	SampleGraphletCumulativeDist(V, Varray, G, k); // This one is faster but less well tested and less well understood.
#endif
	//SampleGraphletUnbiasedMaybe(V, G, k);
	// We should probably figure out a faster sort?
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
#if PARANOID // no point in freeing this stuff since we're about to exit; it can take significant time for large graphs.
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
    if(argc != 4 && argc != 5) Fatal("USAGE: blant [-numCores] {k} {nSamples} {graph.el}. Graph can only have integers as node labels");
    if(argc == 5) // only way to get 5 args is if the user specified "-cores" as the first argument.
    {
	CORES = atoi(argv[1]);
	if(CORES >= 0)
	    Fatal("You specified 5 arguments but the first argument must be specified as -nCores\n"
		"USAGE: blant [-numCores] {k} {nSamples} {graph.el}. Graph can only have integers as node labels");
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
