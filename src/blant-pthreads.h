// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#ifndef BLANT_PTHREADS_H
#define BLANT_PTHREADS_H
#include "blant-fundamentals.h"
#include "sets.h"
#include "graph.h"
#include "tree.h"
#include <stddef.h>

// Define cache line size (typical value for most modern processors)
#define CACHE_LINE_SIZE 64

#if DYNAMIC_CANON_MAP
// Pruning thresholds for the dynamic canonical-map graphlet list. "For now"
// these are fixed: retain the 2,000,000 most frequent graphlets while sampling,
// and report only the top 1,000,000 in the output. Graphlets are ranked by
// count (descending); exact ties are broken by a hash of the canonical ID, so
// the choice among equal-frequency graphlets is not biased by ID magnitude.
#ifndef BLANT_DYNAMIC_PRUNE_KEEP
#define BLANT_DYNAMIC_PRUNE_KEEP 2000000UL
#endif
#ifndef BLANT_DYNAMIC_REPORT_MAX
#define BLANT_DYNAMIC_REPORT_MAX 1000000UL
#endif

// Keep exactly the `keep` most frequent graphlets in *tree and discard the
// rest. This is used twice: during sampling to bound memory, and again at
// output to reduce the survivors to the report size.
void PruneGraphletTree(TREETYPE **tree, size_t keep);
#endif


typedef struct {
    #if !DYNAMIC_CANON_MAP
    double *graphletCount;
    double *graphletConcentration;
    double **graphletDegreeVector;
    double **orbitDegreeVector;
    SET*** communityNeighbors;
    double *canonNumStarMotifs;
    unsigned long *batchRawCount;
    unsigned long batchRawTotalSamples;
    #else
    TREETYPE *graphletCount;
    // Total weight of every sample accumulated into graphletCount, including
    // graphlets that have since been pruned. Used as the concentration
    // denominator so pruning does not distort -F c.
    double graphletTotalWeight;
    #endif
} Accumulators;

// ThreadData structure (now after all required types are defined)
// https://docs.oracle.com/cd/E19120-01/open.solaris/816-5137/tlib-4/index.html
typedef struct ThreadData {
    int k;
    GRAPH *G;
    int varraySize;
    int threadId; // thread number, starting from 0
    long seed;
    Accumulators *accums;

    // batching
    unsigned batchSize;
    unsigned long totalSamples;   // GLOBAL total across all threads
    int totalBatches;             // ceil(totalSamples / batchSize)
} __attribute__((aligned(CACHE_LINE_SIZE))) ThreadData;

// Function declarations
Accumulators* InitializeAccumulatorStruct(GRAPH* G);
void FreeAccumulatorStruct(Accumulators *accums);
void SampleNGraphletsInThreads(long base_seed, int k, GRAPH *G, int varraySize, unsigned numSamples, int numThreads);

// Anytime a function must take an Accumulator as a parameter, but you don't intend on actually using the data, pass it this. 
// The data here is never used, and maintaining only one copy of this saves memory.
extern Accumulators _trashAccumulator;

#endif // BLANT_PTHREADS_H
