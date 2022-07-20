#ifndef CDIJKSTRA_SEEDING_HEADER
#define CDIJKSTRA_SEEDING_HEADER

#include "graph.h"

// Implements a seed line data structure, which maps k-value to the appropriate node. seed_from_file() skips to the given line in the file and reads the seed from that line based on the given graph.

typedef struct {
    unsigned int kval;
    unsigned int node;
} seed_t;

void seed_from_file(seed_t* seed, char* fname, unsigned int line_number, GRAPH* g);

#endif
