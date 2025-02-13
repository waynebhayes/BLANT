#ifndef BLANT_OUTPUT_H
#define BLANT_OUTPUT_H

#include "blant.h"

// ProcessGraphlet returns true if the graphlet was processed, and false if it was removed due to being a duplicate
Boolean ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], const int k, TINY_GRAPH *g, double weight, Accumulators *accums);
Boolean NodeSetSeenRecently(GRAPH *G, unsigned Varray[], int k);

// None of the "Print" routines here actually print anything; they put the string into a constant internal buff, and then
// YOU must print it (immediately is best, since the buffer will be over-written the next time you call the function).
char *PrintGraphletID(Gint_type Gint); // returns a string that YOU must print
char *PrintOrdinal(Gordinal_type GintOrdinal); // returns a string that YOU must print

// For _outputMode = communityDetection, process the graphlet/orbit counts and neighbors
void ProcessNodeOrbitNeighbors(GRAPH *G, Gint_type Gint, Gordinal_type GintOrdinal, unsigned Varray[], TINY_GRAPH *g, int k, double w, unsigned char* perm);
void ProcessNodeGraphletNeighbors(GRAPH *G, Gint_type Gint, Gordinal_type GintOrdinal, unsigned Varray[], TINY_GRAPH *g, int k, double w, unsigned char* perm);

// You should print a node EXCLUSIVELY with this function; it automatically determines if we're supporting names or not.
// If c is non-zero, the character is prepended, using putchar(), before the node is printed (eg ' ', ';' or '\t').
char *PrintNode(char c, int v); // returns a pointer to a constant buffer; YOU need to print it.
// Print a pair of nodes, but sorted so that u<=v, even if they are strings.
char *PrintNodePairSorted(int u, char c, int v); // returns a pointer to a constant buffer, YOU must print it.

char *PrintIndexEntry(Gint_type Gint, Gordinal_type GintOrdinal, unsigned Varray[], int k, double w, unsigned char* perm);
char *PrintIndexOrbitsEntry(Gint_type Gint, Gordinal_type GintOrdinal, unsigned Varray[], int k, double w, unsigned char* perm);

#endif
