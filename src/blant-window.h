#ifndef BLANT_WINDOW_H
#define BLANT_WINDOW_H
#include "blant.h"
#include "heap.h"
#include "graph.h"
#include "multisets.h"
#include <stdbool.h>

//windowRep global Variables
#define WINDOW_SAMPLE_MIN 1 // Find the k-graphlet with minimal canonicalInt
#define WINDOW_SAMPLE_MAX 2 // Find the k-graphlet with maximal canonicalInt
#define WINDOW_SAMPLE_MIN_D 3   // Find the k-graphlet with minimal canonicalInt and balanced numEdges
#define WINDOW_SAMPLE_MAX_D 4   // Find the k-graphlet with maximal canonicalInt and balanced numEdges
#define WINDOW_SAMPLE_LEAST_FREQ_MIN 5 // Find the k-graphlet with least fequent cacnonicalInt. IF there is a tie, pick the minimal one
#define WINDOW_SAMPLE_LEAST_FREQ_MAX 6 // Find the k-graphlet with least fequent cacnonicalInt. IF there is a tie, pick the maximial one
#define WINDOW_SAMPLE_DEG_MAX 7
extern int _windowSampleMethod;

#define WINDOW_LIMIT_UNDEF 0
#define WINDOW_LIMIT_DEGREE 1
#define WINDOW_LIMIT_EDGES 2
extern int _windowRep_limit_method;
extern HEAP * _windowRep_limit_heap;

#define WINDOW_ITER_COMB 1 // using Combination method to sample k-graphlets in Window
#define WINDOW_ITER_DFS 2 // Default way: using DFS-like way.
extern int _windowIterationMethod;

extern int **_windowReps;
extern int _MAXnumWindowRep;
extern int _numWindowRep;
extern int _numWindowRepLimit;
extern int _numWindowRepArrSize;

extern int _topThousandth;
extern int _orbitNumber;
extern char* _odvFile;
extern bool _alphabeticTieBreaking;

extern int _windowSize;
extern Boolean _window;
extern SET *_windowRep_allowed_ambig_set;
extern int _windowRep_min_num_edge;
extern float *_graphNodeImportance;
extern Boolean _supportNodeImportance;

extern Boolean _windowRep_limit_neglect_trivial;

void FindWindowRepInWindow(GRAPH *G, SET *W, int *windowRepInt, int *D, char perm[]);
void ProcessWindowRep(GRAPH *G, int *VArray, int windowRepInt);
void ProcessWindowDistribution(GRAPH *G, SET *V, unsigned Varray[], int k, TINY_GRAPH *prev_graph, SET *prev_node_set, SET *intersect_node);

#endif
