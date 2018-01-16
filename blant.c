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
static TINY_GRAPH *TinyGraphInducedFromGraph(GRAPH *G, SET *V)
{
    unsigned array[MAX_TSET], nV = SetToArray(array, V), i, j;
    assert(nV <= MAX_TSET);
    TINY_GRAPH *Gv = TinyGraphAlloc(nV);
    for(i=0; i < nV; i++) for(j=i+1; j < nV; j++)
        if(GraphAreConnected(G, array[i], array[j]))
            TinyGraphConnect(Gv, i, j);
    return Gv;
}

// Given the big graph G and an integer k, return a k-graphlet from G.
// Caller is responsible for freeing the SET that is returned.
static SET *SampleGraphlet(GRAPH *G, int k)
{
    int edge = G->numEdges * drand48(), i, v1, v2;
    SET *V = SetAlloc(G->n), *outSet = SetAlloc(G->n);
    int nOut = 0, outbound[G->n]; // vertices one step outside the boundary of V
    v1 = G->edgeList[2*edge];
    v2 = G->edgeList[2*edge+1];
    SetAdd(V, v1);
    SetAdd(V, v2);
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
		; // must terminate since k < G->n
	    outbound[nOut++] = j;
	    j = 0;
	}
	else
	    j = nOut * drand48();
	v1 = outbound[j];
	SetDelete(outSet, v1);
	SetAdd(V, v1);
	outbound[j] = outbound[--nOut];	// nuke v1 from the list of outbound by moving the last one to its place
	for(j=0; j<G->degree[v1];j++)
	{
	    v2 = G->neighbor[v1][j];
	    if(!SetIn(outSet, v2) && !SetIn(V, v2))
		SetAdd(outSet, (outbound[nOut++] = v2));
	}
    }
    assert(SetCardinality(V) == k);
    SetFree(outSet);
    return V;
}

static SET *SampleGraphletUnbiased(GRAPH *G, int k)
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
    assert(SetCardinality(V) == k);
    // Now check that it's connected
    GRAPH *graphette = GraphInduced(NULL, G, V);
    if(GraphBFS(graphette, 0, k, nodeArray, distArray) < k)
    {
	GraphFree(graphette);
	SetFree(V);
	return SampleGraphletUnbiased(G, k);
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

// This is the single-core version of blant. It used to be the main() function, which is why it takes (argc,argv[]).
// Now it simply gets called once for each core we use when run with multiple cores.
int blant(int argc, char *argv[])
{
    int k=atoi(argv[1]), i, j;
    _k = k;
    int numSamples = atoi(argv[2]);
    FILE *fp = fopen(argv[3], "r");
    GRAPH *G = GraphReadEdgeList(fp, true); // sparse = true
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

    for(i=0; i<numSamples; i++)
    {
	SET *V = SampleGraphlet(G, k);
	//SET *V = SampleGraphletUnbiased(G, k);
	TINY_GRAPH *g = TinyGraphInducedFromGraph(G, V);
	int Gint = TinyGraph2Int(g,k);
	for(j=0;j<k;j++) perm[j]=0;
	ExtractPerm(perm, Gint);
	unsigned Varray[maxK+1];
	SetToArray(Varray, V);
	//printf("K[%d]=%d [%d];", Gint, K[Gint], canon_list[K[Gint]]);
	printf("%d", K[Gint]); // Note this is the ordinal of the canonical, not its bit representation
	for(j=0;j<k;j++) printf(" %d", Varray[(unsigned)perm[j]]);
	puts("");
	TinyGraphFree(g);
	SetFree(V);
    }
    GraphFree(G);
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
