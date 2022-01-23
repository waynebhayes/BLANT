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
    {"g1file", required_argument, NULL, '1'},
    {"g2file", required_argument, NULL, '2'},
    {"delta", required_argument, NULL, 'd'},
    {0, 0, 0, 0},
};

int main(int argc, char** argv) {
    string graph1_fname;
    string graph2_fname;
    double delta;

    // TODO: switch to long options where appropriate
    char option;
    int option_index = -1;

    while ((option = getopt_long(argc, argv, "d:", long_opts, &option_index)) != -1) {
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

        case '?':
            exit(ERROR_INVALID_OPTION);

        default:
            exit(ERROR_UNKNOWN);
        }
    }

    cout << "graph 1 filename: " << graph1_fname << endl;
    cout << "graph 2 filename: " << graph2_fname << endl;
    cout << "delta: " << delta << endl;

    return 0;
}
