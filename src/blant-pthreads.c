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
#include "tree.h"

atomic_u64_t nextIndex CACHE_ALIGNED = 0; // single definition

// Worker thread stack size: 1 MiB
// Empirical measurements show ~257 bytes actual
// 1 MiB provides ~4000x safety margin for larger k values and different networks
static const size_t kWorkerThreadStackBytes = 1 * 1024 * 1024;

#if DYNAMIC_CANON_MAP
// The graphlet-frequency map is keyed by a canonical graphlet ID (Gint_type),
// which is 128 bits when TINY_SET_SIZE == 16, so it does not fit directly in a
// foint. Store keys as heap pointers and give the tree a dereferencing
// comparator plus deep copy/free functions (same pattern as copyInt128Ptr and
// freeInt128Ptr in blant-utils.c).
static int CmpGintKey(foint a, foint b) {
    const Gint_type x = *(const Gint_type *)a.v;
    const Gint_type y = *(const Gint_type *)b.v;
    return (x < y) ? -1 : (x > y) ? 1 : 0;
}
static foint CopyGintKey(foint v) {
    Gint_type *p = Malloc(sizeof(Gint_type));
    *p = *(const Gint_type *)v.v;
    return (foint){.v = p};
}
static void FreeGintKey(foint v) { free(v.v); }

// Starting value for the tiebreak hash. Overriding it (at compile time)
// changes which graphlets win an exact frequency tie, but not any count.
#ifndef BLANT_DYNAMIC_TIEBREAK_SEED
#define BLANT_DYNAMIC_TIEBREAK_SEED 0xcbf29ce484222325ULL // FNV-1a offset basis
#endif

// Hash a canonical graphlet ID to 64 bits. Only used to break exact frequency
// ties. The canonical ID is a structural encoding (it is the minimum adjacency
// matrix over all relabelings), so ordering by it would systematically favour
// graphlets that happen to encode to small numbers; hashing removes that bias.
// Byte-wise FNV-1a plus a splitmix64 finalizer works for 32-, 64- and
// 128-bit Gint_type without knowing which one is in use.
static uint64_t _gintHash(Gint_type x) {
    const unsigned char *p = (const unsigned char *)&x;
    uint64_t h = BLANT_DYNAMIC_TIEBREAK_SEED;
    size_t i;
    for (i = 0; i < sizeof x; ++i) {
        h ^= p[i];
        h *= 0x100000001b3ULL; // FNV-1a prime
    }
    h ^= h >> 30;
    h *= 0xbf58476d1ce4e5b9ULL;
    h ^= h >> 27;
    h *= 0x94d049bb133111ebULL;
    h ^= h >> 31;
    return h;
}

// One accumulated graphlet and its (weighted) sample count. File-local: only
// the pruning helper needs to rank graphlets, so this type stays out of the
// headers.
typedef struct {
    Gint_type gint;
    double count;
    // Hash of gint, used to break frequency ties without the structural bias of
    // the numeric canonical identifier (which is what gint itself would be).
    uint64_t tiebreak;
} GraphletCountEntry;

// Gather state for _collectCounts.
// List, index pair
typedef struct {
    GraphletCountEntry *entries;
    size_t n;
} CollectState;

//helper to modify entries in TreeTraverse
static int _collectGraphletCount(foint globals, foint key, foint data) {
    CollectState *st = (CollectState *)globals.v;
    GraphletCountEntry *e = &st->entries[st->n++];
    e->gint = *(Gint_type *)key.v;
    e->count = data.d;
    e->tiebreak = _gintHash(e->gint);
    return 1;
}

// Descending count, then ascending hash of the canonical ID, then ascending
// canonical ID as a last-resort collision fallback. The hash breaks exact
// frequency ties without the structural bias of the raw canonical ID.
static int _countDescCmp(const void *pa, const void *pb) {
    const GraphletCountEntry *a = (const GraphletCountEntry *)pa;
    const GraphletCountEntry *b = (const GraphletCountEntry *)pb;
    if (a->count > b->count) return -1;
    if (a->count < b->count) return 1;
    if (a->tiebreak < b->tiebreak) return -1;
    if (a->tiebreak > b->tiebreak) return 1;
    if (a->gint < b->gint) return -1;
    if (a->gint > b->gint) return 1;
    return 0;
}

// Copy `tree` into a malloc'd array ordered by descending count (see
// _countDescCmp). The caller must Free() *out; *out is NULL for an empty tree.
static void _collectCounts(TREETYPE *tree, GraphletCountEntry **out) {
    *out = NULL;
    if (!tree || tree->n == 0) return;
    // tree->n is the number of live distinct keys, so one exact allocation is
    // enough; a duplicate key does not increment it.
    size_t n = tree->n;
    GraphletCountEntry *entries = Malloc(n * sizeof(*entries));
    CollectState st = {entries, 0};
    TreeTraverse((foint){.v = &st}, tree, _collectGraphletCount);
    if (st.n != n)
        Fatal("graphlet tree count mismatch: tree->n = %zu, traversed %zu", n,
              st.n);
    qsort(entries, n, sizeof(*entries), _countDescCmp);
    *out = entries;
}

void PruneGraphletTree(TREETYPE **tree, size_t keep) {
    if (!tree || !*tree || keep == 0) return;
    if ((*tree)->n <= keep) return;
    GraphletCountEntry *entries = NULL;
    _collectCounts(*tree, &entries); // count equals (*tree)->n (checked there)
    TREETYPE *pruned = TreeAlloc(CmpGintKey, CopyGintKey, FreeGintKey, NULL, NULL);
    size_t i;
    for (i = 0; i < keep; ++i) {
        Gint_type gint = entries[i].gint;
        foint key = {.v = &gint}, info = {.d = entries[i].count};
        TreeInsert(pruned, key, info);
    }
    TreeFree(*tree);
    *tree = pruned;
    Free(entries);
}
#endif

Accumulators* InitializeAccumulatorStruct(GRAPH* G) {
    // Ensure size is a multiple of alignment
    size_t size = sizeof(Accumulators);
    size_t aligned_size = ((size + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE) * CACHE_LINE_SIZE;
    
    Accumulators *accums = aligned_alloc(CACHE_LINE_SIZE, aligned_size);

    if (!accums) {
        perror("Failed to allocate memory for Accumulators");
        exit(EXIT_FAILURE);
    }

    memset(accums, 0, aligned_size); // zero out everything including padding
    #if !DYNAMIC_CANON_MAP
    accums->graphletCount = Calloc(_numCanon, sizeof(double));
    accums->graphletConcentration = Calloc(_numCanon, sizeof(double));
    accums->batchRawCount = Calloc(_numCanon, sizeof(unsigned long));
    accums->batchRawTotalSamples = 0;
    accums->orbitDegreeVector = Calloc(_numOrbits, sizeof(double *));
    accums->canonNumStarMotifs = Malloc(_numCanon * sizeof(double));
    for (int i = 0; i < _numCanon; i++) accums->canonNumStarMotifs[i] = -1;
    
    // initialize GDV vectors if needed
    if(_outputMode & outputGDV || (_outputMode & communityDetection && _communityMode=='g')) {
        accums->graphletDegreeVector = Malloc(_numCanon * sizeof(double*));
        if (!accums->graphletDegreeVector) Fatal("Failed to allocate GDV memory");
        for(int i = 0; i < _numCanon; i++) {
            size_t bytes = G->n * sizeof(double);
            // Ensure that size is a multiple of CACHE_LINE_SIZE
            if (bytes % CACHE_LINE_SIZE != 0) {
                bytes = (bytes / CACHE_LINE_SIZE + 1) * CACHE_LINE_SIZE;
            }
            accums->graphletDegreeVector[i] = aligned_alloc(CACHE_LINE_SIZE, bytes);
            if (!accums->graphletDegreeVector[i]) {
                perror("Failed to allocate memory for GDV vector");
                exit(EXIT_FAILURE);
            }
            memset(accums->graphletDegreeVector[i], 0, bytes);
        }
    }

    // initialize ODV vectors if needed
    if(_outputMode & outputODV || (_outputMode & communityDetection && _communityMode=='o')) {
        for(int i=0; i<_numOrbits; i++) {
            size_t bytes = G->n * sizeof(double);
            // Ensure that size is a multiple of CACHE_LINE_SIZE
            if (bytes % CACHE_LINE_SIZE != 0) {
                bytes = (bytes / CACHE_LINE_SIZE + 1) * CACHE_LINE_SIZE;
            }
            accums->orbitDegreeVector[i] = aligned_alloc(CACHE_LINE_SIZE, bytes);
            if (!accums->orbitDegreeVector[i]) {
                perror("Failed to allocate memory for ODV vector");
                exit(EXIT_FAILURE);
            }
            memset(accums->orbitDegreeVector[i], 0, bytes);
        }
    }
    if(_outputMode & communityDetection) accums->communityNeighbors = (SET***) calloc(G->n, sizeof(SET**));
    #else
    accums->graphletCount = TreeAlloc(CmpGintKey, CopyGintKey, FreeGintKey, NULL, NULL);
    #endif
    return accums;
}

void FreeAccumulatorStruct(Accumulators *accums) {
    #if !DYNAMIC_CANON_MAP
    Free(accums->graphletCount);
    Free(accums->graphletConcentration);
    Free(accums->canonNumStarMotifs);
    Free(accums->batchRawCount);
    if(_outputMode & outputGDV || (_outputMode & communityDetection && _communityMode=='g')) {
        for (int i=0; i<_numCanon; i++) free(accums->graphletDegreeVector[i]);
        Free(accums->graphletDegreeVector);
    }
    if(_outputMode & outputODV || (_outputMode & communityDetection && _communityMode=='o')) {
        for(int i=0; i<_numOrbits; i++) free(accums->orbitDegreeVector[i]);
        Free(accums->orbitDegreeVector);
    }
    if(_outputMode & communityDetection) free(accums->communityNeighbors);
    #else
    if(accums->graphletCount) TreeFree(accums->graphletCount);
    #endif
    free(accums);
}
#if !DYNAMIC_CANON_MAP
void SampleNGraphletsInThreads(long base_seed, int k, GRAPH *G, int varraySize, unsigned numSamples, int numThreads) {
    // Reset the global “next” counter before launching workers
    ATOMIC_STORE_U64(&nextIndex, 0);

    if (numThreads < 1) numThreads = 1;
    if (numSamples < 0) numSamples = 0;

    pthread_t threads[numThreads];
    ThreadData threadData[numThreads];


    // Choose a batch size
    unsigned batchSize = G->numEdges * sqrt(G->n);
    if (batchSize <= 0) batchSize = 1;

    if (numSamples > 0 && batchSize > numSamples) batchSize = numSamples;

    // make sure no single thread claims the entire workload when numSamples is small
    int fairShare = numSamples > 0 ? (numSamples + numThreads - 1) / numThreads : 1;
    if (batchSize > fairShare) batchSize = fairShare;
    if (batchSize <= 0) batchSize = 1;

    int totalBatches = batchSize > 0 ? (numSamples + batchSize - 1) / batchSize : 1; // Ceiling division to cover all samples
    if (totalBatches <= 0) totalBatches = 1;

    // Setup pthread attributes for custom stack size
    pthread_attr_t attr;
    pthread_attr_t *attr_ptr = NULL;
    Boolean use_custom_stack = false;
    
    if (pthread_attr_init(&attr) == 0) {
        if (pthread_attr_setstacksize(&attr, kWorkerThreadStackBytes) == 0) {
            use_custom_stack = true;
            attr_ptr = &attr;
        } else {
            Warning("Failed to set worker thread stack size to %zu bytes; using system default.", kWorkerThreadStackBytes);
            pthread_attr_destroy(&attr);
        }
    } else {
        Warning("Failed to initialize pthread attributes; using system default stack size.");
    }

    Note("SampleN base_seed = %ld", base_seed);
    // initialize the threads and their data
    for (int t = 0; t < numThreads; t++)
    {
        threadData[t].k = k;
        threadData[t].G = G;
        threadData[t].varraySize = varraySize;
        threadData[t].threadId = t;
        threadData[t].seed = base_seed + t;
        threadData[t].accums = NULL; // Will be initialized in the thread

        // batching params consumed by RunBlantInThread
        threadData[t].batchSize    = batchSize;
        threadData[t].totalSamples = (unsigned long)numSamples;
        threadData[t].totalBatches = totalBatches;

        // create thread with error handling
        if (pthread_create(&threads[t], attr_ptr, RunBlantInThread, &threadData[t]) != 0) {
            Fatal("Failed to create thread");
        }
    }
    
    // Clean up pthread attributes if we used them
    if (use_custom_stack) {
        pthread_attr_destroy(&attr);
    }

    // wait for each thread to finish execution, then accumulate data from the thread into the passed accumulator
    for (unsigned t = 0; t < numThreads; t++) {
        pthread_join(threads[t], NULL);
        for (int i = 0; i < _numCanon; i++) {
            _graphletConcentration[i] += threadData[t].accums->graphletConcentration[i];
            _graphletCount[i] += threadData[t].accums->graphletCount[i];
            if (_canonNumStarMotifs[i] == -1) _canonNumStarMotifs[i] = threadData[t].accums->canonNumStarMotifs[i];
            // Accumulate batch counters
            _batchRawCount[i] += threadData[t].accums->batchRawCount[i];
        }
        // Accumulate total batch samples
        _batchRawTotalSamples += threadData[t].accums->batchRawTotalSamples;

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
#endif