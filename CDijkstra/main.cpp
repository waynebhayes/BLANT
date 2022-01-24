#include <getopt.h>
#include <stdlib.h>

#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

#define ERROR_INVALID_OPTION 1
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

    string graph1_seed_line;
    string graph2_seed_line;

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
            graph1_seed_line = optarg;
            break;

        case 'M':
            graph2_seed_line = optarg;
            break;

        case '?':
            exit(ERROR_INVALID_OPTION);

        default:
            exit(ERROR_UNKNOWN);
        }
    }

    // TODO: enforce required options

    cout << "graph 1 file: " << graph1_fname << endl;
    cout << "graph 2 file: " << graph2_fname << endl;
    cout << "delta: " << delta << endl;
    cout << "sim file: " << sim_fname << endl;
    cout << "runs: " << runs << endl;

    return 0;
}
