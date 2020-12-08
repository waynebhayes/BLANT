#ifndef BLANT_OUTPUT_H
#define BLANT_OUTPUT_H

#include "blant.h"

// ProcessGraphlet returns true if the graphlet was processed, and false if it was removed due to being a duplicate
Boolean ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], const int k, char perm[], TINY_GRAPH *g);
void PrintCanonical(int GintOrdinal);
Boolean NodeSetSeenRecently(GRAPH *G, unsigned Varray[], int k);

// You should print a node EXCLUSIVELY with this function; it automatically determines if we're supporting names or not.
// If c is non-zero, the character is appended, using putchar(), after the node is printed (usually '\t', or ' ' or '\n');
void PrintNode(int v, char c);

#endif
