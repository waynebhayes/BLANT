#include <getopt.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <utility>

#include "graph.h"
#include "matrix.h"
#include "seeding.h"
#include "sim.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

using cdijkstra::get_seed_line;
using cdijkstra::get_sim;
using cdijkstra::Graph;
using cdijkstra::Matrix;

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
    string graph1_fname;
    string graph2_fname;

    string sim_fname;

    string graph1_seed_fname;
    string graph2_seed_fname;

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
    if (graph1_fname.empty()) {
        cerr << "error: graph1 file is required" << endl;
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (graph2_fname.empty()) {
        cerr << "error: graph2 file is required" << endl;
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (sim_fname.empty()) {
        cerr << "error: similarity file is required" << endl;
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

/*
    if (graph1_seed_fname.empty()) {
        cerr << "error: graph1 seed file is required" << endl;
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (graph2_seed_fname.empty()) {
        cerr << "error: graph2 seed file is required" << endl;
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (graph1_seed_line == -1) {
        cerr << "error: graph1 seed line is required" << endl;
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }

    if (graph2_seed_line == -1) {
        cerr << "error: graph2 seed line is required" << endl;
        exit(ERROR_INVALID_OR_MISSING_OPTION);
    }
*/
    Graph g1(graph1_fname);
    Graph g2(graph2_fname);

    Matrix<double> sim = get_sim(sim_fname, g1, g2);

    // pair<unsigned int, unsigned int> g1_seed_line = get_seed_line(graph1_seed_fname, graph1_seed_line, g1);
    // pair<unsigned int, unsigned int> g2_seed_line = get_seed_line(graph2_seed_fname, graph2_seed_line, g2);
    /*
    // check kvals
    if (g1_seed_line.first != g2_seed_line.first) {
        cerr << "error: seed k-values don't match" << endl;
        exit(ERROR_UNKNOWN);
    }
    */

    return 0;
}
