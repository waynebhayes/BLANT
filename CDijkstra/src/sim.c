#include "sim.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void sim_set(sim_t* sim, int r, int c, double val) {
    int i = sim_index(sim, r, c);
    sim->mat[i] = val;
}

void sim_from_file(sim_t* sim, char* fname, GRAPH* g1, GRAPH* g2) {
    // each line looks like:
    // NODE1 NODE2 SIM
    FILE* sim_file = fopen(fname, "r");
    if (sim_file == NULL) {
        return;
    }

    double* mat = malloc(sizeof(double) * g1->n * g2->n);
    if (mat == NULL) {
        return;
    }

    sim->mat = mat;
    sim->rows = g1->n;
    sim->cols = g2->n;

    char* buf = NULL;
    ssize_t bufsize = 0;

    while (getline(&buf, &bufsize, sim_file) != -1) {
        char* first_delim = strchr(buf, ' ');
        *first_delim = '\0';
        int n1 = GraphNodeName2Int(g1, buf);

        char* second_delim = strchr(first_delim + 1, ' ');
        *second_delim = '\0';
        int n2 = GraphNodeName2Int(g2, first_delim + 1);

        double sim_val = atof(second_delim + 1);
        sim_set(sim, n1, n2, sim_val);
    }
}

void sim_free(sim_t* sim) {
    free(sim->mat);
}

int sim_index(sim_t* sim, int i, int j) {
    return i * sim->rows + j;
}

double sim_get(sim_t* sim, int r, int c) {
    int i = sim_index(sim, r, c);
    return sim->mat[i];
}
