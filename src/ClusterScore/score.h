#define NDEBUG 1
#define PARANOID_ASSERTS 0
#include "misc.h"
#include "graph.h"
#include "sets.h"
#include "rand48.h"
#include <math.h>
#include <signal.h>

typedef struct _clustering {
    GRAPH *G;
    SET **clusters;
    unsigned nC, *clusSize;
} CLUSTERING;

// Allocate a clustering and return it with the zero'th cluster assigned to all nodes in G.
CLUSTERING *ClusteringAlloc(GRAPH *G);
unsigned ClusterEdgeCount(GRAPH *G, SET *c, unsigned num);
double ScoreOneCluster(GRAPH *G, SET *c, unsigned n);
double ScoreClustering(CLUSTERING *C);
void SplitCluster(CLUSTERING *C, unsigned who);
Boolean TrySplit(CLUSTERING *C);
Boolean TryMerge(CLUSTERING *C);
Boolean TryMove(CLUSTERING *C);
void OutputClustering(int sig);
Boolean ScoreIteration();
