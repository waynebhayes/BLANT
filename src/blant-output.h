#ifndef BLANT_OUTPUT_H
#define BLANT_OUTPUT_H

#include "blant.h"

void ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], const int k, char perm[], TINY_GRAPH *g);
void PrintCanonical(int GintOrdinal);
Boolean GraphletSeenRecently(GRAPH *G, unsigned Varray[], int k);

#endif