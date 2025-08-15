#include "blant.h"
#include "blant-output.h"
#include "blant-sampling.h"
#include "blant-pthreads.h"
#include <pthread.h>


// it would be potentially more efficient to just use a static set of accumulators struct across each batch, then memset(0) them at the start of every batch
Accumulators* InitializeAccumulatorStruct(GRAPH* G) {
    // memory NEEDS to be optimized- but to do this we need to understand Ocalloc. Why is it used here? What's its purpose? Why not use malloc?
    Accumulators *accums = calloc(1, sizeof(Accumulators)); // I don't really know how libwayne's Omalloc and Ofree works, confer during meeting later
    // initialize GDV vectors if needed
    if(_outputMode & outputGDV || (_outputMode & communityDetection && _communityMode=='g'))
        for(int i=0; i<_numCanon; i++) accums->graphletDegreeVector[i] = calloc(G->n, sizeof(**accums->graphletDegreeVector));
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
    pthread_t threads[numThreads];
    ThreadData threadData[numThreads];
    int samplesPerThread = numSamples / numThreads;
    int leftover = numSamples - (samplesPerThread * numThreads);
    // seed the threads with a base seed that may or may not be specified
    long base_seed = seed == -1 ? GetFancySeed(false) : seed;

    // initialize the threads and their data
    for (unsigned t = 0; t < numThreads; t++)
    {
        threadData[t].samplesPerThread = samplesPerThread;
        // distribute the leftover samples
        if (leftover-- > 0) threadData[t].samplesPerThread++;
        threadData[t].k = k;
        threadData[t].G = G;
        threadData[t].varraySize = varraySize;
        threadData[t].threadId = t;
        threadData[t].seed = base_seed + t; // each thread has it's own unique seed
        threadData[t].accums = InitializeAccumulatorStruct(G);

        pthread_create(&threads[t], NULL, RunBlantInThread, &threadData[t]);
    }

    // wait for each thread to finishe execution, then accumulate data from the thead into the passed accumulator
    for (unsigned t = 0; t < _numThreads; t++) {
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
