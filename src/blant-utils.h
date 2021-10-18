#ifndef BLANT_UTILS_H
#define BLANT_UTILS_H

#include "blant.h"

typedef struct nwd { // this struct is needed for enumerateNodeOrder
    int node;
    int heur;
} node_wheur;

Boolean arrayIn(int* arr, int size, int item);
void ExtractPerm(char perm[_k], int i);
void InvertPerm(char inverse[_k], const char perm[_k]);
TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, int *Varray);
int getMaximumIntNumber(int K);
int asccompFunc(const foint i, const foint j);
int descompFunc(const void *a, const void *b);
int nwh_descompFunc(const void *a, const void *b);
int nwh_asccompFunc(const void *a, const void *b);
void SetGlobalCanonMaps(void);
void LoadMagicTable();
int* enumerateNodeOrder(GRAPH *G);
void enumerateNodeOrderHelper(GRAPH *G, node_wheur* orderArray, int start, int end, int layer);
void antidupFillNextStepSet(SET **next_step, SET **deg_set, GRAPH *G, int *prev_nodes_array, int prev_nodes_count, int *node_order);

#endif
