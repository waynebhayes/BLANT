#ifndef IMPORTANCE_H
#define IMPORTANCE_H


#include "graph.h"
#include "blant-utils.h" // to get node_wdegree

#define IMPORTANCE_DEG 10

typedef struct {
    int **lists;
    int *sizes;
} ADJ_LIST;

void getDoubleDegreeArr(double *double_degree_arr, GRAPH *G);
int *enumerateImportanceDegreeOrder(GRAPH *G);
void getImportances(double *importances, GRAPH *G);
ADJ_LIST *generateAdjacencyList(GRAPH *G);
void freeAdjacencyList(ADJ_LIST *adjList, GRAPH *G);
void fillSortedNodes(int *sortedNodes, GRAPH *G);
void removeFromAdjacencyList(ADJ_LIST *adjList, int base, int neigh);


#endif
