// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#ifndef BLANT_PTHREADS_H
#define BLANT_PTHREADS_H
#include "blant-fundamentals.h"
#include "sets.h"
#include "graph.h"


typedef struct {
    // local accumulator values, they function the same as globals but ARE LOCAL TO THREADS
    double graphletCount[MAX_CANONICALS];
    double graphletConcentration[MAX_CANONICALS];
    double *graphletDegreeVector[MAX_CANONICALS];
    double *orbitDegreeVector[MAX_ORBITS];
    SET*** communityNeighbors;
    double canonNumStarMotifs[MAX_CANONICALS];
} Accumulators;

Accumulators* InitializeAccumulatorStruct(GRAPH* G);
void FreeAccumulatorStruct(Accumulators *accums);

// Anytime a function must take an Accumulator as a parameter, but you don't intend on actually using the data, pass it this. 
// The data here is never used, and maintaining only one copy of this saves memory.
extern Accumulators _trashAccumulator;

// https://docs.oracle.com/cd/E19120-01/open.solaris/816-5137/tlib-4/index.html
typedef struct {
    int samplesPerThread;
    int k;
    GRAPH *G;
    int varraySize;
    int threadId; // thread number, starting from 0
    long seed;
    Accumulators *accums;
} ThreadData;

void SampleNGraphletsInThreads(int seed, int k, GRAPH *G, int varraySize, int numSamples, int numThreads);


#endif // BLANT_PTHREADS_H