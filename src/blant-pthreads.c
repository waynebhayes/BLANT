// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include "blant.h"
#include "blant-output.h"
#include "blant-sampling.h"
#include "blant-pthreads.h"
#include "rand48.h"
#include "misc.h"
#include <pthread.h>
#include <errno.h>
#include "atomic_utils.h"
atomic_u64_t nextIndex CACHE_ALIGNED = 0; // single definition

Accumulators* InitializeAccumulatorStruct(GRAPH* G) {

    Accumulators *accums = aligned_alloc(CACHE_LINE_SIZE, sizeof(Accumulators));

    if (!accums) {
        perror("Failed to allocate memory for Accumulators");
        exit(EXIT_FAILURE);
    }

    memset(accums, 0, sizeof(Accumulators)); // zero out everything
    
    // initialize GDV vectors if needed
    if(_outputMode & outputGDV || (_outputMode & communityDetection && _communityMode=='g')) {
        accums->graphletDegreeVector = malloc(_numCanon * sizeof(double*));
        if (!accums->graphletDegreeVector) Fatal("Failed to allocate GDV memory");
        for(int i = 0; i < _numCanon; i++) {
            accums->graphletDegreeVector[i] = calloc(G->n, sizeof(double));
            if (!accums->graphletDegreeVector[i]) Fatal("Failed to allocate GDV entry memory");
        }
    }

    // initialize ODV vectors if needed
    if(_outputMode & outputODV || (_outputMode & communityDetection && _communityMode=='o'))
        for(int i=0; i<_numOrbits; i++) accums->orbitDegreeVector[i] = calloc(G->n, sizeof(**accums->orbitDegreeVector));
    // initialize communityNeighbors if needed
    if(_outputMode & communityDetection) accums->communityNeighbors = (SET***) calloc(G->n, sizeof(SET**));
        
    for (int i = 0; i < _numCanon; i++) accums->canonNumStarMotifs[i] = -1; // initialize to -1, meaning "not yet initialized"

    return accums;
}

void FreeAccumulatorStruct(Accumulators *accums) {
    assert(_numCanon <= MAX_CANONICALS);
    if(_outputMode & outputGDV || (_outputMode & communityDetection && _communityMode=='g'))
        for (int i=0; i<_numCanon; i++) free(accums->graphletDegreeVector[i]);
    if(_outputMode & outputODV || (_outputMode & communityDetection && _communityMode=='o'))
        for(int i=0; i<_numOrbits; i++) if (accums->orbitDegreeVector[i] != NULL) free(accums->orbitDegreeVector[i]);
    if(_outputMode & communityDetection) free(accums->communityNeighbors);
    free(accums);
}

void SampleNGraphletsInThreads(int seed, int k, GRAPH *G, int varraySize, int numSamples, int numThreads) {
    // Reset the global “next” counter before launching workers
    ATOMIC_STORE_U64(&nextIndex, 0);

    if (numThreads < 1) numThreads = 1;
    if (numSamples < 0) numSamples = 0;

    pthread_t threads[numThreads];
    ThreadData threadData[numThreads];

    // Choose a batch size
    int batchSize = G->numEdges * sqrt(G->n) * sqrt(_numThreads);
    if (batchSize > numSamples) batchSize = numSamples > 0 ? numSamples : 1;

    int totalBatches = (numSamples + batchSize - 1) / batchSize; // Ceiling division to cover all samples
    if (totalBatches <= 0) totalBatches = 1;

    // seed the threads with a base seed that may or may not be specified
    long base_seed = seed == -1 ? GetFancySeed(false) : seed;

    // initialize the threads and their data
    for (int t = 0; t < numThreads; t++)
    {
        threadData[t].k = k;
        threadData[t].G = G;
        threadData[t].varraySize = varraySize;
        threadData[t].threadId = t;
        threadData[t].seed = base_seed + t;
        threadData[t].accums = InitializeAccumulatorStruct(G);

        // batching params consumed by RunBlantInThread
        threadData[t].batchSize    = batchSize;
        threadData[t].totalSamples = (unsigned long)numSamples;
        threadData[t].totalBatches = totalBatches;

        // create thread with error handling
        if (pthread_create(&threads[t], NULL, RunBlantInThread, &threadData[t]) != 0) {
            Fatal("Failed to create thread");
        }
    }

    // wait for each thread to finish execution, then accumulate data from the thread into the passed accumulator
    for (unsigned t = 0; t < numThreads; t++) {
        pthread_join(threads[t], NULL);
        for (int i = 0; i < _numCanon; i++) {
            _graphletConcentration[i] += threadData[t].accums->graphletConcentration[i];
            _graphletCount[i] += threadData[t].accums->graphletCount[i];
            if (_canonNumStarMotifs[i] == -1) _canonNumStarMotifs[i] = threadData[t].accums->canonNumStarMotifs[i];
        }

        if (_outputMode & outputODV || (_outputMode & communityDetection && _communityMode=='o')) {
            for(int i=0; i<_numOrbits; i++) {
            for(int j=0; j<G->n; j++) {
                _orbitDegreeVector[i][j] += threadData[t].accums->orbitDegreeVector[i][j];
            }
            }
        }
        if (_outputMode & outputGDV || (_outputMode & communityDetection && _communityMode=='g')) {
            for(int i=0; i<_numCanon; i++) {
            for(int j=0; j<G->n; j++) {
                _graphletDegreeVector[i][j] += threadData[t].accums->graphletDegreeVector[i][j];
            }
            }
        }
        if (_outputMode & communityDetection) {
            int numCommunities = (_communityMode=='o') ? _numOrbits : _numCanon;
            for(int i=0; i<G->n; i++) {
                if(threadData[t].accums->communityNeighbors[i]) {
                if(!_communityNeighbors[i]) {
                    _communityNeighbors[i] = (SET**) Calloc(numCommunities, sizeof(SET*));
                }
                for(int j=0; j<numCommunities; j++) {
                    if(threadData[t].accums->communityNeighbors[i][j]) {
                    if(!_communityNeighbors[i][j]) {
                        _communityNeighbors[i][j] = SetAlloc(G->n);
                    }
                    _communityNeighbors[i][j] = SetUnion(_communityNeighbors[i][j], _communityNeighbors[i][j], threadData[t].accums->communityNeighbors[i][j]);
                    }
                }
                }
            }
        }
    }

    // In each threadData[t] accumulator, the GDV and ODV vectors are allocated with Ocalloc, and must be freed with Ofree
    // However, if anything else has been allocated with Ocalloc or Omalloc with libwayne BETWEEN the time this function starts and ends
    // same goodbye to that memory and say hello to segfault -Ethan
    for (unsigned t = 0; t < numThreads; t++) {
        FreeAccumulatorStruct(threadData[t].accums);
    }
}
