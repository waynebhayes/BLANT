#ifndef BLANT_H
#define BLANT_H

#include "tinygraph.h"
#include "sets.h"
#include "blant-window.h"
#include "graph.h"

#define USE_MarsenneTwister 0
#if USE_MarsenneTwister
#include "libwayne/MT19937/mt19937.h"
static MT19937 *_mt19937;
#define RandomSeed SeedMt19937
void SeedMt19937(int seed) {
    if(_mt19937) Mt19937Free(_mt19937);
    _mt19937 = Mt19937Alloc(seed);
}
double RandomUniform(void) {
    return Mt19937NextDouble(_mt19937);
}
#else
#include "rand48.h"
#define RandomUniform drand48
#define RandomSeed srand48
#endif

#define GEN_SYN_GRAPH 0

// This is the maximum graphlet size that BLANT supports.  Cannot be bigger than 8.
// Currently only used to determine the amount of static memory to allocate.
#define maxK 8

#define maxBk (1 << (maxK*(maxK-1)/2)) // maximum number of entries in the canon_map

#define mcmc_d 2 // arbitrary d graphlet size < k for MCMC algorithm. Should always be 2 or k-1

#define MAX_CANONICALS	12346	// This is the number of canonical graphettes for k=8
extern int _numCanon, _canonList[MAX_CANONICALS];
extern int _numConnectedCanon;

#define MAX_ORBITS	79264	// This is the number of orbits for k=8

// BLANT represents a graphlet using one-half of the adjacency matrix (since we are assuming symmetric, undirected graphs)
// We have a choice of using the upper or lower triangle. We prefer the lower triangle because that's what Jesse uses
// (the graphlet / orbit generation code of Ine Melckenbeeck and friends Ghent university), and they published first.
// NOTE THAT IF YOU CHOOSE UPPER TRIANGLE THEN THE TESTS IN THE MAKEFILE WILL FAIL.
#define LOWER_TRIANGLE	1

// Once we find which canonical graphlet corresponds to a sampled graphlet, we want to know the permutation between the
// two.  We default to the permutation from the canonical to the sampled non-canonical; thus, when we list the nodes
// in the graphlet on BLANT's output, the ordering means that the columns of identical graphlets correspond to an exact
// local alignment between the nodes listed.  Listing them in the other direction doesn't seem to have much use.
#define PERMS_CAN2NON	1

#define CANON_DIR "canon_maps"
//#define CANON_DIR "/var/preserve/Graphette/canon_maps" // if you happen to put it there...

#define DEFAULT_BLANT_DIR "."
extern char* _BLANT_DIR;

#define PARANOID_ASSERTS 1	// turn on paranoid checking --- slows down execution by a factor of 2-3

extern unsigned long int *_graphletDegreeVector[MAX_CANONICALS];
extern unsigned long int    *_orbitDegreeVector[MAX_ORBITS];
extern double *_doubleOrbitDegreeVector[MAX_ORBITS];

// If you're squeemish then use this one to access the degrees:
#define ODV(node,orbit)       _orbitDegreeVector[orbit][node]
#define GDV(node,graphlet) _graphletDegreeVector[graphlet][node]

// Enable the code that uses C++ to parse input files?
#define SHAWN_AND_ZICAN 0
extern unsigned int _Bk;

extern int _numOrbits, _orbitList[MAX_CANONICALS][maxK];
extern int _alphaList[MAX_CANONICALS];

extern char **_nodeNames, _supportNodeNames;
extern unsigned int _k;
extern short int *_K;
extern SET *_connectedCanonicals;

enum OutputMode {undef, indexGraphlets, indexOrbits, indexMotifs, indexMotifOrbits,
    kovacsPairs, kovacsAllOrbits, graphletFrequency, outputODV, outputGDV,
    graphletDistribution // used in Windowing
};
extern enum OutputMode _outputMode;
extern int _outputMapping[MAX_CANONICALS];

extern unsigned long int _graphletCount[MAX_CANONICALS];
extern int **_graphletDistributionTable;
extern double _graphletConcentration[MAX_CANONICALS];

enum CanonicalDisplayMode {undefined, ordinal, decimal, binary, orca, jesse};
extern enum CanonicalDisplayMode _displayMode;
extern int _orbitCanonMapping[MAX_ORBITS]; // Maps orbits to canonical (including disconnected)

enum FrequencyDisplayMode {freq_display_mode_undef, count, concentration};
extern enum FrequencyDisplayMode _freqDisplayMode;

extern int _orca_orbit_mapping[58]; // Mapping from orbit indices in orca_ordering to our orbits
extern int _connectedOrbits[MAX_ORBITS];
extern int _numConnectedOrbits;

extern unsigned int _numConnectedComponents;
extern unsigned int *_componentSize; // number of nodes in each CC
extern unsigned int *_whichComponent; // will be an array of size G->n specifying which CC each node is in.
extern SET **_componentSet;

extern double *_cumulativeProb;


#define SPARSE true // do not try false at the moment, it's broken

Boolean NodeSetSeenRecently(GRAPH *G, unsigned Varray[], int k);
void SampleGraphlet(GRAPH *G, SET *V, unsigned Varray[], int k);

void BuildGraph(TINY_GRAPH* G, int Gint);
int TinyGraph2Int(TINY_GRAPH *g, int numNodes);
short int* mapCanonMap(char* BUF, short int *K, int k);
SET *canonListPopulate(char *BUF, int *canon_list, int k); // returns a SET containing list of connected ordinals
int orbitListPopulate(char *BUF, int orbit_list[MAX_CANONICALS][maxK],  int orbit_canon_mapping[MAX_ORBITS], int numCanon, int k);
void orcaOrbitMappingPopulate(char *BUF, int orca_orbit_mapping[58], int k);
char** convertToEL(char* file); // from convert.cpp
void BuildGraphletsAsSeeds(GRAPH* G, SET* prev_nodes, int numSamplesPerNode, int *tempCountPtr);

#endif
