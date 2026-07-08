// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#ifndef __SCORE_H
#define __SCORE_H 1
//#define NDEBUG 1
//#define PARANOID_ASSERTS 0
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
Boolean ScoreIteration(double pBad, Boolean running);
double ScoreClustering(CLUSTERING *C);
void SplitCluster(CLUSTERING *C, unsigned who);
Boolean TrySplit(CLUSTERING *C);
Boolean TryMerge(CLUSTERING *C);
Boolean TryMove(CLUSTERING *C);
void OutputClustering(int sig);

CLUSTERING *_C;
double _fullScore;
int _largestCluster;
#endif
