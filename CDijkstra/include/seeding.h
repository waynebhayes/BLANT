#ifndef CDIJKSTRA_SEEDING_HEADER
#define CDIJKSTRA_SEEDING_HEADER

#include "graph.h"

typedef struct {
    unsigned int kval;
    unsigned int node;
} seed_t;

void seed_from_file(seed_t* seed, char* fname, unsigned int line_number, GRAPH* g);

#endif
