#ifndef BLANT_SAMPLING_H
#define BLANT_SAMPLING_H
#include "blant.h"
#include "heap.h"
#include "graph.h"
#include "queue.h"
#include "multisets.h"

#define SAMPLE_FAYE 6
#define MAX_TRIES 100		// max # of tries in cumulative sampling before giving up
#define ALLOW_DISCONNECTED_GRAPHLETS 0
#define PARANOID_ASSERTS 1	// turn on paranoid checking --- slows down execution by a factor of 2-3

#ifndef RESERVOIR_MULTIPLIER
// this*k is the number of steps in the Reservoir walk. 8 seems to work best, empirically.
#define RESERVOIR_MULTIPLIER 8
#endif

// Below are the sampling methods
#define SAMPLE_FROM_FILE 0
#define SAMPLE_ACCEPT_REJECT 1	// makes things REALLY REALLY slow.  Like 10-100 samples per second rather than a million.
#define SAMPLE_NODE_EXPANSION 2	// sample using uniform node expansion; about 100,000 samples per second
#define SAMPLE_EDGE_EXPANSION 3	// Fastest, up to a million samples per second
#define SAMPLE_RESERVOIR 4	// Lu Bressan's reservoir sampler, reasonably but not entirely unbiased.
#define SAMPLE_MCMC 5 // MCMC Algorithm estimates graphlet frequency with a random walk

extern int _sampleMethod;
extern FILE *_sampleFile; // if _sampleMethod is SAMPLE_FROM_FILE
extern char _sampleFileEOF;
extern unsigned long int _acceptRejectTotalTries;
extern int _samplesPerEdge;
extern int _numSamples;

extern unsigned _MCMC_L; // walk length for MCMC algorithm. k-d+1 with d almost always being 2.
extern Boolean _MCMC_EVERY_EDGE; // Should MCMC restart at each edge

static SET *SampleGraphletNodeBasedExpansion(SET *V, int *Varray, GRAPH *G, int k, int whichCC);
static SET *SampleGraphletFaye(SET *V, int *Varray, GRAPH *G, int k, int whichCC);
static SET *SampleGraphletFromFile(SET *V, int *Varray, GRAPH *G, int k);
static SET *SampleGraphletEdgeBasedExpansion(SET *V, int *Varray, GRAPH *G, int k, int whichCC);
static SET *SampleGraphletLuBressanReservoir(SET *V, int *Varray, GRAPH *G, int k, int whichCC);
static SET *SampleGraphletAcceptReject(SET *V, int *Varray, GRAPH *G, int k);
static SET *SampleGraphletMCMC(SET *V, int *Varray, GRAPH *G, int k, int whichCC);
static SET *SampleGraphletLuBressan_MCMC_MHS_without_Ooze(SET *V, int *Varray, GRAPH *G, int k);
static SET *SampleGraphletLuBressan_MCMC_MHS_with_Ooze(SET *V, int *Varray, GRAPH *G, int k);
static int NumReachableNodes(TINY_GRAPH *g, int startingNode);
void WalkLSteps(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G, int k, int cc, int edge);
void SampleGraphlet(GRAPH *G, SET *V, unsigned Varray[], int k);
void initializeSlidingWindow(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G, int windowSize, int edge);
void crawlOneStep(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G);
int *MCMCGetNeighbor(int *Xcurrent, GRAPH *G);

#endif
