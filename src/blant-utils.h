#ifndef BLANT_UTILS_H
#define BLANT_UTILS_H

#include "blant.h"

typedef struct nwd { // this struct is needed for enumerateDegreeOrder
    int node;
    int degree_val;
} node_wdegree;

void ExtractPerm(char perm[_k], int i);
TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, int *Varray);
int getMaximumIntNumber(int K);
int asccompFunc(const foint i, const foint j);
int descompFunc(const void *a, const void *b);
int nwd_descompFunc(const void *a, const void *b);
void SetGlobalCanonMaps(void);
void LoadMagicTable();
int* enumerateDegreeOrder(GRAPH *G);
void enumerateDegreeOrderHelper(GRAPH *G, node_wdegree* orderArray, int start, int end, int layer);

#endif
