#ifndef BLANT_SAMPLING_H
#define BLANT_SAMPLING_H
#include "blant.h"
#include "heap.h"
#include "graph.h"
#include "queue.h"
#include "multisets.h"
#include "blant-utils.h"

#define MAX_TRIES 100		// max # of tries in cumulative sampling before giving up
#define ALLOW_DISCONNECTED_GRAPHLETS 0

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
#define SAMPLE_FAYE 6
#define SAMPLE_INDEX 7 // Use deterministic walk to find seeds which are used for extensions
#define SAMPLE_MCMC_EC 8 // cover all edges in G at least once, then stop
#define SAMPLE_KRMCMC 9 // MCMC Algorithm but it restarts each time it gets a graphlet

extern int _sampleMethod, _sampleSubmethod;
extern FILE *_sampleFile; // if _sampleMethod is SAMPLE_FROM_FILE
extern char _sampleFileEOF;
extern unsigned long int _acceptRejectTotalTries;
extern int _samplesPerEdge;
extern double _g_overcount;

extern unsigned _MCMC_L; // walk length for MCMC algorithm. k-d+1 with d almost always being 2.
extern Boolean _MCMC_EVERY_EDGE; // Should MCMC restart at each edge
extern GRAPH *_EDGE_COVER_G;

// Each of these samples a graphlet and return its weight (default 1 unless _weighted)
double SampleGraphletNodeBasedExpansion(GRAPH *G, SET *V, unsigned *Varray, int k, int whichCC);
double SampleGraphletFaye(GRAPH *G, SET *V, unsigned *Varray, int k, int whichCC);
double SampleGraphletFromFile(GRAPH *G, SET *V, unsigned *Varray, int k);
double SampleGraphletEdgeBasedExpansion(GRAPH *G, SET *V, unsigned *Varray, int k, int whichCC);
double SampleGraphletLuBressanReservoir(GRAPH *G, SET *V, unsigned *Varray, int k, int whichCC);
double SampleGraphletAcceptReject(GRAPH *G, SET *V, unsigned *Varray, int k);
double SampleGraphletMCMC(GRAPH *G, SET *V, unsigned *Varray, int k, int whichCC);
double SampleGraphletLuBressan_MCMC_MHS_without_Ooze(GRAPH *G, SET *V, unsigned *Varray, int k);
double SampleGraphletLuBressan_MCMC_MHS_with_Ooze(GRAPH *G, SET *V, unsigned *Varray, int k);
void SampleGraphletIndexAndPrint(GRAPH* G, unsigned *prev_nodes_array, int prev_nodes_count, double *heur_arr); // returns void instead of double because this function isn't called in SampleGraphlet
void WalkLSteps(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G, int k, int cc, int edge);
double SampleGraphlet(GRAPH *G, SET *V, unsigned Varray[], int k, int cc); // call with cc=G->n to allow unbiased choice
void initializeSlidingWindow(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G, int windowSize, int edge);
void crawlOneStep(MULTISET *XLS, QUEUE *XLQ, int* X, GRAPH *G);
int *MCMCGetNeighbor(int *Xcurrent, GRAPH *G);

#endif
