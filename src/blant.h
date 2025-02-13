#ifndef BLANT_H
#define BLANT_H

#include "blant-fundamentals.h" // defining k related constants, including SELF_LOOPS

#include "misc.h"
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
extern double _confidence; // for confidence intervals
extern Boolean _earlyAbort;  // Can be set true by anybody anywhere, and they're responsible for producing a warning as to why

#define mcmc_d 2 // arbitrary d graphlet size < k for MCMC algorithm. Should always be 2 or k-1

extern unsigned long _known_canonical_count[]; //known number of canonicals for k=0 to 12 inclusive (13 entries)
	//{0, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168, 1018997864, 165091172592};
	//  k=1  2  3   4   5   6     7     8      9        10         11          12

// This ugly code is in preparation for allowing k>8 lookup tables (using associative arrays)
#if TINY_SET_SIZE == 16
  #if __GLIBCXX_TYPE_INT_N_0 == 128
    typedef __uint128 Gint_type; // need 120 bits for undirected adjacency matrix for k=16...
    #define GINT_FMT "%ull"
    typedef __uint128 Gordinal_type; // ... and max numCanonicals = 64001015704527557894928 < 2^76... what type is that???
    #define GORDINAL_FMT "%ull"
    #define MAX_BINTREE_K 16
  #elif long_width >= 120
    typedef unsigned long Gint_type, Gordinal_type; //... numCanonicals = 64001015704527557894928 < 2^76
    #define GINT_FMT "%lu"
    #define GORDINAL_FMT "%ul"
    #define MAX_BINTREE_K 16
  #elif long_width >= 55
    typedef unsigned long Gint_type; // for k=11, need 55 bits in undirected adjacency matrix, so 64 bits sufficient...
    #define GINT_FMT "%lu"
    #if int_width >=30
      typedef unsigned Gordinal_type; // ... and max numCanonicals = 1018997864 < 2^30, so 32 bits are sufficient
      #define GORDINAL_FMT "%u"
    #else
      typedef unsigned long Gordinal_type;
      #define GORDINAL_FMT "%ul"
    #endif
    #define MAX_BINTREE_K 11
  #else
    #error "cannot do TINY_SET_SIZE 16 due to no integers being wide enough"
  #endif
#elif TINY_SET_SIZE == 8
  #if int_width >= 28 && short_width >= 16
    typedef unsigned Gint_type; // at k=8, max lookup index is 2^28, so we need 32 bits...
    #define GINT_FMT "%u"
    typedef unsigned short Gordinal_type; //... and max numCanonicals is 12348 < 2^16, so 16 bits is sufficient
    #define GORDINAL_FMT "%hu"
    #define MAX_BINTREE_K 8
  #else
    #error "cannot do TINY_SET_SIZE 8 due to no integers long enough"
  #endif
#else
    #error "unknwon TINY_SET_SIZE"
#endif

extern Gordinal_type _numCanon, _numConnectedCanon;
extern char _canonNumEdges[MAX_CANONICALS];
extern double _totalStarMotifs;
extern Gint_type _canonList[MAX_CANONICALS];

double GetCPUseconds(void);
void Int2TinyGraph(TINY_GRAPH* G, Gint_type Gint);
Gint_type TinyGraph2Int(TINY_GRAPH *g, int numNodes);
Gordinal_type * mapCanonMap(char* BUF, Gordinal_type *K, int k);
SET *canonListPopulate(char *BUF, Gint_type *canon_list, int k, char *canon_num_edges); // returns a SET containing list of connected ordinals
Gint_type orbitListPopulate(char *BUF, Gint_type orbit_list[MAX_CANONICALS][MAX_K],
    Gordinal_type orbit_canon_mapping[MAX_ORBITS], char orbit_canon_node_mapping[MAX_ORBITS], Gordinal_type numCanon, int k);
void orcaOrbitMappingPopulate(char *BUF, int orca_orbit_mapping[58], int k);
char** convertToEL(char* file); // from convert.cpp

#define CANON_DIR "canon_maps"
//#define CANON_DIR "/var/preserve/Graphette/canon_maps" // if you happen to put it there...

#define DEFAULT_BLANT_DIR "."
extern const char* _BLANT_DIR;

extern double *_graphletDegreeVector[MAX_CANONICALS];
extern double       *_orbitDegreeVector[MAX_ORBITS];
extern double *_doubleGraphletDegreeVector[MAX_CANONICALS];
extern double *_doubleOrbitDegreeVector[MAX_ORBITS], _absoluteCountMultiplier;

// If you're squeemish then use this one to access the degrees:
#define ODV(node,orbit)       _orbitDegreeVector[orbit][node]
#define GDV(node,graphlet) _graphletDegreeVector[graphlet][node]

// Enable the code that uses C++ to parse input files?
#define SHAWN_AND_ZICAN 0
extern unsigned int _Bk;

extern Gint_type _numOrbits, _orbitList[MAX_CANONICALS][MAX_K], _alphaList[MAX_CANONICALS];

extern char **_nodeNames, _supportNodeNames;
extern unsigned int _k, _min_edge_count;
extern Gordinal_type *_K; // works because max numCanonicals = 12348 < 2^16, but will need to be > 16 bits for k>8.
extern Gordinal_type L_K_Func(Gint_type Gint);
#if DYNAMIC_CANON_MAP
#define L_K(Gint) L_K_Func(Gint)
#else
#define L_K(Gint) (_K ? _K[Gint] : L_K_Func(Gint))
#endif
extern SET *_connectedCanonicals;

// When in community mode: for each node, and for each orbit it is in, a set of its neigbors when it apppears at that orbit
extern SET ***_communityNeighbors;
extern char _communityMode; // 'g' for graphlet or 'o' for orbit; default ???

// Use POWERS OF 2 so we can combine them.
enum OutputMode {undef=0, indexGraphlets=1, indexGraphletsRNO=2, indexOrbits=4, indexMotifs=8, communityDetection=16,
    indexMotifOrbits=32, predict=64, predict_merge=128, graphletFrequency=256, outputODV=512, outputGDV=1024,
    graphletDistribution=2048 // used in Windowing
};
extern enum OutputMode _outputMode; // note they can be LOGINAL OR'd together; modes can overlap!
extern int _outputMapping[MAX_CANONICALS], _canonNumStarMotifs[MAX_CANONICALS];

extern double _graphletCount[MAX_CANONICALS];
extern int **_graphletDistributionTable;
extern double _graphletConcentration[MAX_CANONICALS];
extern unsigned long _batchRawCount[MAX_CANONICALS], _batchRawTotalSamples; // batches for confidence intervals
enum PrecisionMode {PrecUndef, mean, worst};
enum PrecisionWeights {PrecWtNone, PrecWtRaw, PrecWtLog};

enum CanonicalDisplayMode {undefined, ordinal, decimal, binary, orca, jesse, noncanonical};
extern enum CanonicalDisplayMode _displayMode;
extern Gordinal_type _orbitCanonMapping[MAX_ORBITS]; // Maps orbits to canonical (ordinal value, including disconnected graphlets)
extern char _orbitCanonNodeMapping[MAX_ORBITS]; // Maps orbits to canonical nodes

enum FrequencyDisplayMode {freq_display_mode_undef, freq_display_mode_count, freq_display_mode_concentration, freq_display_mode_estimate_absolute};
extern enum FrequencyDisplayMode _freqDisplayMode;

extern int _orca_orbit_mapping[58]; // Mapping from orbit indices in orca_ordering to our orbits
extern int _connectedOrbits[MAX_ORBITS];
extern int _numConnectedOrbits;

extern int _numConnectedComponents;
extern int *_componentSize; // number of nodes in each CC
extern int *_whichComponent; // will be an array of size G->n specifying which CC each node is in.
extern SET **_componentSet;

// from -i "nodes of interest" file, the list of nodes to start all graphlet samples, or the identity mapping if not
extern int *_startNodes, _numStartNodes;
extern SET *_startNodeSet;

extern double *_cumulativeProb;
extern Boolean _child, _weighted;
extern Boolean _rawCounts;
extern int _quiet; // suppress notes and/or warnings, higher value = more quiet

#define SPARSE true // do not try false at the moment, it's broken

Boolean NodeSetSeenRecently(GRAPH *G, unsigned Varray[], int k);


// Headers for multithreading
typedef struct {
    // local accumulator values, they function the same as globals but ARE LOCAL TO THREADS
    double graphletCount[MAX_CANONICALS];
    double graphletConcentration[MAX_CANONICALS];
} Accumulators;

// https://docs.oracle.com/cd/E19120-01/open.solaris/816-5137/tlib-4/index.html
typedef struct {
    int numSamples;
    int k;
    GRAPH *G;
    int varraySize;
    Accumulators accums;
} ThreadData;
void* RunBlantInThread(void* arg);

#endif
