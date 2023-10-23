#ifndef BLANT_H
#define BLANT_H

#include "blant-fundamentals.h" // defining k related constants, including SELF_LOOPS

#include "sets.h"
#include "tinygraph.h"
#include "blant-window.h"
#include "graph.h"
#include "Oalloc.h"
// #include "mem-debug.h" // need this if you use ENABLE_MEMORY_TRACKING() at the top of main().

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

#define USE_INSERTION_SORT 0
#define GEN_SYN_GRAPH 0

#define MAX_POSSIBLE_THREADS 64 // set this to something reasonable on your machine (eg odin.ics.uci.edu has 64 cores)
extern int _JOBS, _MAX_THREADS;
extern unsigned long _numSamples;
extern Boolean _earlyAbort;  // Can be set true by anybody anywhere, and they're responsible for producing a warning as to why

#define mcmc_d 2 // arbitrary d graphlet size < k for MCMC algorithm. Should always be 2 or k-1

#define MAX_CANONICALS	12346	// This is the number of canonical graphettes for k=8

// This ugly code is in preparation for allowing k>8 lookup tables (using associative arrays)
#if TINY_SET_SIZE == 16
#if __GLIBCXX_TYPE_INT_N_0 == 128
    typedef __uint128 Gint_type;
    #define MAX_BINTREE_K 16
#else
    typedef unsigned long long Gint_type;
    #define MAX_BINTREE_K 11
#endif
#elif TINY_SET_SIZE == 8
typedef unsigned Gint_type;
    #define MAX_BINTREE_K 8
#else
#error "unknwon TINY_SET_SIZE"
#endif

extern int _numCanon, _numConnectedCanon, _canonNumEdges[MAX_CANONICALS];
extern Gint_type _canonList[MAX_CANONICALS];

#define MAX_ORBITS	79264	// This is the number of orbits for k=8

void Int2TinyGraph(TINY_GRAPH* G, Gint_type Gint);
Gint_type TinyGraph2Int(TINY_GRAPH *g, int numNodes);
short int* mapCanonMap(char* BUF, short int *K, int k);
SET *canonListPopulate(char *BUF, Gint_type *canon_list, int k, int *canon_num_edges); // returns a SET containing list of connected ordinals
int orbitListPopulate(char *BUF, int orbit_list[MAX_CANONICALS][MAX_K],  int orbit_canon_mapping[MAX_ORBITS],
    int orbit_canon_node_mapping[MAX_ORBITS], int numCanon, int k);
void orcaOrbitMappingPopulate(char *BUF, int orca_orbit_mapping[58], int k);
char** convertToEL(char* file); // from convert.cpp

#define CANON_DIR "canon_maps"
//#define CANON_DIR "/var/preserve/Graphette/canon_maps" // if you happen to put it there...

#define DEFAULT_BLANT_DIR "."
extern const char* _BLANT_DIR;

#define PARANOID_ASSERTS 0	// turn on copious assert checking --- slows down execution by a factor of 2-3

extern double *_graphletDegreeVector[MAX_CANONICALS];
extern double    *_orbitDegreeVector[MAX_ORBITS];
extern double *_doubleOrbitDegreeVector[MAX_ORBITS], _absoluteCliqueCount, _absoluteCountMultiplier;

// If you're squeemish then use this one to access the degrees:
#define ODV(node,orbit)       _orbitDegreeVector[orbit][node]
#define GDV(node,graphlet) _graphletDegreeVector[graphlet][node]

// Enable the code that uses C++ to parse input files?
#define SHAWN_AND_ZICAN 0
extern unsigned int _Bk;

extern int _numOrbits, _orbitList[MAX_CANONICALS][MAX_K];
extern int _alphaList[MAX_CANONICALS];

extern char **_nodeNames, _supportNodeNames;
extern unsigned int _k, _min_edge_count;
extern short int *_K;
extern unsigned int L_K_Func(Gint_type Gint);
#define L_K(Gint) (_K ? _K[Gint] : L_K_Func(Gint))
extern SET *_connectedCanonicals;

// When in community mode: for each node, and for each orbit it is in, a set of its neigbors when it apppears at that orbit
extern SET ***_communityNeighbors;
extern char _communityMode; // 'g' for graphlet or 'o' for orbit; default ???

enum OutputMode {undef, indexGraphlets, indexGraphletsRNO, indexOrbits, indexMotifs, communityDetection,
    indexMotifOrbits, predict, predict_merge, graphletFrequency, outputODV, outputGDV,
    graphletDistribution // used in Windowing
};
extern enum OutputMode _outputMode;
extern int _outputMapping[MAX_CANONICALS];

extern double _graphletCount[MAX_CANONICALS];
extern int **_graphletDistributionTable;
extern double _graphletConcentration[MAX_CANONICALS];

enum CanonicalDisplayMode {undefined, ordinal, decimal, binary, orca, jesse};
extern enum CanonicalDisplayMode _displayMode;
extern int _orbitCanonMapping[MAX_ORBITS]; // Maps orbits to canonical (ordinal value, including disconnected graphlets)
extern int _orbitCanonNodeMapping[MAX_ORBITS]; // Maps orbits to canonical nodes

enum FrequencyDisplayMode {freq_display_mode_undef, count, concentration};
extern enum FrequencyDisplayMode _freqDisplayMode;

extern int _orca_orbit_mapping[58]; // Mapping from orbit indices in orca_ordering to our orbits
extern int _connectedOrbits[MAX_ORBITS];
extern int _numConnectedOrbits;

extern int _numConnectedComponents;
extern int *_componentSize; // number of nodes in each CC
extern int *_whichComponent; // will be an array of size G->n specifying which CC each node is in.
extern SET **_componentSet;

extern double *_cumulativeProb;
extern Boolean _child, _weighted;

#define SPARSE true // do not try false at the moment, it's broken

Boolean NodeSetSeenRecently(GRAPH *G, unsigned Varray[], int k);

#endif
