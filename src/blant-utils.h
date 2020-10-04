#ifndef BLANT_UTILS_H
#define BLANT_UTILS_H

#include "../blant.h"

void ExtractPerm(char perm[_k], int i);
TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, int *Varray);
int getMaximumIntNumber(int K);
int asccompFunc(const foint i, const foint j);
int descompFunc(const void *a, const void *b);
void SetGlobalCanonMaps(void);
void LoadMagicTable();

#endif