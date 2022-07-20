#ifndef CDIJKSTRA_SIM_HEADER
#define CDIJKSTRA_SIM_HEADER

#include <string.h>
#include "graph.h"

// A node similarity support data structure, basically a n1 x n2 matrix of similarity values. Implemented as a single buffer in row-major order, because that's faster than allocating for each row.

typedef struct {
    int rows;
    int cols;
    double* mat; // the contents of the matrix, stored in row-major order
} sim_t;

void sim_from_file(sim_t* sim, char* fname, GRAPH* g1, GRAPH* g2);
void sim_free(sim_t* sim);

double sim_get(sim_t* sim, int r, int c);

#endif
