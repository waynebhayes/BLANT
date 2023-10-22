#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

// #include "graph.h"
// #include "matrix.h"
// #include "seeding.h"
// #include "sim.h"

#include "graph.h"
#include "seeding.h"
#include "sim.h"

#define ERROR_INVALID_OR_MISSING_OPTION 1
#define ERROR_UNKNOWN -1

static struct option long_opts[] = {
    {"graph1", required_argument, NULL, '1'},
    {"graph2", required_argument, NULL, '2'},
    {"delta", required_argument, NULL, 'd'},
    {"sim", required_argument, NULL, 's'},
    {"runs", required_argument, NULL, 'r'},
    {"graph1-seeds", required_argument, NULL, 'S'},
    {"graph2-seeds", required_argument, NULL, 'T'},
    {"graph1-seedline", required_argument, NULL, 'L'},
    {"graph2-seedline", required_argument, NULL, 'M'},
    {0, 0, 0, 0},
};

int main(int argc, char** argv) {
    char* graph1_fname = "";
    char* graph2_fname = "";

    char* sim_fname = "";

    char* graph1_seed_fname = "";
    char* graph2_seed_fname = "";

    unsigned int graph1_seed_line = -1;
    unsigned int graph2_seed_line = -1;

    double delta;
    int runs = 1;

    char option;
    int option_index = -1;

    while ((option = getopt_long(argc, argv, "d:r:s:", long_opts, &option_index)) != -1) {
        switch (option) {
        case '1':
            graph1_fname = optarg;
            break;

        case '2':
            graph2_fname = optarg;
            break;

        case 'd':
            delta = atof(optarg);
            break;

        case 's':
            sim_fname = optarg;
            break;

        case 'r':
            runs = atoi(optarg);
            break;

        case 'S':
            graph1_seed_fname = optarg;
            break;

        case 'T':
            graph2_seed_fname = optarg;
            break;

        case 'L':
            graph1_seed_line = atoi(optarg);
            break;

        case 'M':
            graph2_seed_line = atoi(optarg);
            break;

        case '?':
            exit(ERROR_INVALID_OR_MISSING_OPTION);

        default:
            exit(ERROR_UNKNOWN);
        }
    }

    // TODO: enforce required options
    if (strlen(graph1_fname) == 0) {
        fprintf(stderr, "error: graph1 file is required\n");
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (strlen(graph2_fname) == 0) {
        fprintf(stderr, "error: graph2 file is required\n");
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (strlen(sim_fname) == 0) {
        fprintf(stderr, "error: similarity file is required\n");
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (strlen(graph1_seed_fname) == 0) {
        fprintf(stderr, "error: graph1 seed file is required\n");
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (strlen(graph2_seed_fname) == 0) {
        fprintf(stderr, "error: graph2 seed file is required\n");
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (graph1_seed_line == -1) {
        fprintf(stderr, "error: graph1 seed line is required\n");
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (graph2_seed_line == -1) {
        fprintf(stderr, "error: graph2 seed line is required\n");
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    FILE* g1_file = fopen(graph1_fname, "r");
    FILE* g2_file = fopen(graph2_fname, "r");

    GRAPH* g1 = GraphReadEdgeList(g1_file, 1, 1,false);
    GRAPH* g2 = GraphReadEdgeList(g2_file, 1, 1,false);

    // NOTE: just for testing
    GraphPrintConnections(stdout, g1);
    printf("\n");
    GraphPrintConnections(stdout, g2);

    sim_t sim;
    sim_from_file(&sim, sim_fname, g1, g2);

    int i, j;

    for (i = 0; i < sim.rows; i++) {
        char* n1_name = g1->name[i];

        for (j = 0; j < sim.cols; j++) {
            char* n2_name = g2->name[j];
            double sim_val = sim_get(&sim, i, j);
            
            if (sim_val > 0) {
                printf("%s %s %f\n", n1_name, n2_name, sim_val);
            }
        }
    }

    seed_t g1_seed, g2_seed;
    seed_from_file(&g1_seed, graph1_seed_fname, graph1_seed_line, g1);
    seed_from_file(&g2_seed, graph2_seed_fname, graph2_seed_line, g2);

    if (g1_seed.kval != g2_seed.kval) {
        fprintf(stderr, "error: seed k-values don't match\n");
        exit(ERROR_UNKNOWN);
    }

    return 0;
}
