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
#include "blant.h"
#include "sets.h"

char USAGE[] = "synthetic -k k Gtarget.el Gsynth.el blant.Gt.index blant.Gs.index\n - output is new Gsynth.el\n";

Boolean _supportNodeNames = false;

// The following is the most compact way to store the permutation between a non-canonical and its canonical representative,
// when k=8: there are 8 entries, and each entry is a integer from 0 to 7, which requires 3 bits. 8*3=24 bits total.
// For simplicity we use the same 3 bits per entry, and assume 8 entries, even for k<8.  It wastes memory for k<4, but
// makes the coding much simpler.
typedef unsigned char kperm[3]; // The 24 bits are stored in 3 unsigned chars.

static unsigned int _Bk, _k; // _k is the global variable storing k; _Bk=actual number of entries in the canon_map for given k.
static int _numCanon, _canonList[MAX_CANONICALS];
static int _numOrbits, _orbitList[MAX_CANONICALS][maxK];

// Here's where we're lazy on saving memory, and we could do better.  We're going to allocate a static array
// that is big enough for the 256 million permutations from non-canonicals to canonicals for k=8, even if k<8.
// So we're allocating 256MBx3=768MB even if we need much less.  I figure anything less than 1GB isn't a big deal
// these days. It needs to be aligned to a page boundary since we're going to mmap the binary file into this array.
static kperm Permutations[maxBk] __attribute__ ((aligned (8192)));
// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
static short int _K[maxBk] __attribute__ ((aligned (8192)));


// Assuming the global variable _k is set properly, go read in and/or mmap the big global
// arrays related to canonical mappings and permutations.
void SetGlobalCanonMaps(void)
{
    assert(3 <= _k && _k <= 8);
    _Bk = (1 <<(_k*(_k-1)/2));
    char BUF[BUFSIZ];
    int i;
    _numCanon = canonListPopulate(BUF, _canonList, _k);
    _numOrbits = orbitListPopulate(BUF, _orbitList, _k);
    mapCanonMap(BUF, _K, _k);
    sprintf(BUF, CANON_DIR "/perm_map%d.bin", _k);
    int pfd = open(BUF, 0*O_RDONLY);
    kperm *Pf = Mmap(Permutations, _Bk*sizeof(Permutations[0]), pfd);
    assert(Pf == Permutations);
}


int main(int argc, char *argv[])
{
    int i, opt, numSamples=0, j, line;

    if(argc == 1)
    {
	fprintf(stderr, "%s\n", USAGE);
	exit(1);
    }

    _k = 0;

    while((opt = getopt(argc, argv, "k:")) != -1)
    {
	switch(opt)
	{
	case 'k': _k = atoi(optarg);
	    if(!(3 <= _k && _k <= 8)) Fatal("k must be between 3 and 8\n%s", USAGE);
	    break;
	default: Fatal("unknown option %c\n%s", opt, USAGE);
	}
    }
    SetGlobalCanonMaps(); // needs _k to be set

    GRAPH *G[2]; // G[0] is the target, G[1] is the synthetic.
    for(i=0; i<2; i++)
    {
	if(!argv[optind]) Fatal("no input graph file specified\n%s", USAGE);
	FILE *fpGraph = fopen(argv[optind], "r");
	if(!fpGraph) Fatal("cannot open graph input file '%s'\n", argv[optind]);
	// Read it in using native Graph routine.
	G[i] = GraphReadEdgeList(fpGraph, true, _supportNodeNames); // sparse=true
	if(_supportNodeNames)
	    assert(G[i]->name);
	fclose(fpGraph);
	optind++;
    }

    // The distribution of graphlets (squiggly plot vectors)
    int D[2][_numCanon];
    for(i=0; i<_numCanon; i++) D[0][i] = D[1][i] = 0;

    int nSamples[2], S, **BLANT[2];
    for(i=0;i<2;i++)
    {
	char cmd[BUFSIZ];
	// Get number of samples
	sprintf(cmd, "wc -l < %s", argv[optind]);
	FILE *fp = popen(cmd, "r");
	assert(fp);
	fscanf(fp, "%d", &nSamples[i]);
	assert(nSamples[i] > 0);
	pclose(fp);
	fp=fopen(argv[optind], "r");
	assert(fp);
	BLANT[i] = (int**)Malloc(nSamples[i] * sizeof(int*));
	for(line=0; line<nSamples[i]; line++)
	{
	    BLANT[i][line] = (int*)Malloc((_k+1)*sizeof(int));
	    for(j=0; j<_k+1; j++)
		assert(1 == fscanf(fp, "%d", &(BLANT[i][line][j])));
	    ++D[i][BLANT[i][line][0]]; // squiggly plot
	}
	assert(fscanf(fp, "%d", &line) < 1); // ensure there's nothing left to read.
	fclose(fp);
	optind++;
    }
    assert(nSamples[0] == nSamples[1]);
    S = nSamples[0];
    for(i=0;i<2;i++)
    {
	int testCount = 0;
	for(j=0; j<_numCanon; j++) testCount += D[i][j];
	assert(testCount == S);
    }

    SET **samples = (SET**)Malloc(G[1]->n * sizeof(SET*));
    for(i=0; i < G[1]->n; i++)
	samples[i] = SetAlloc(S);
    for(line = 0; line < S; line++)
	for(j=1; j<=_k; j++)
	    SetAdd(samples[BLANT[1][line][j]], line);

/*

while(not done---either some number of iterations, or objective function says we're too far away)
{
    pick an edge e = (v1, v2) from synthetic network Gs; // we're going to move it elsewhere
    do {eNew = random u1 and u2} while(IsEdge(Gs, u1,u2)); // find a non-edge
    // delete edge e from Gs and add eNew, then update the index + distribution vectors
    assert(GraphIsEdge(Gs, v1, v2));
    GraphDisconnect(Gs, v1, v2); // remove edge e from Gs
    for(L in samples[v1])
        if(SetIn(samples[v2], L)) // both v1 and v2 are in sample L
            decrement counter for Ds[line L's first column]
            Recompute canonical and set the first column appropriately
            increment counter for Ds[line L's new first column]
        end if // we've now determined the new graphlet for line L of Gs's index after removing edge e
    end for
    GraphConnect(Gs, u1, u2); // add edge eNew to Gs
    Now do the same thing for(L in samples[u1]) and if(SetIn(samples[u2], L)....
    ..... bla bla // determined new graphlet for line L of Gs's index after adding edge eNew
    end for
    Compute objective function comparing distribution(Gs) to distribution(Gt)
    // note we should do this "efficiently" using a delta rather than recomputing whole objective
    if(we need to revert the change)
        revert all of the above: move eNew back to e and reverse all the changes.
    end if
}
output edge list of Gs
*/

}
