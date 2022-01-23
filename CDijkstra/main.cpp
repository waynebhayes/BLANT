#include <getopt.h>
#include <stdlib.h>

#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

int main(int argc, char** argv) {
    string graph1_fname;
    string graph2_fname;
    double delta;

    // TODO: switch to long options where appropriate
    char option;

    while ((option = getopt(argc, argv, "p:q:d:")) != -1) {
        switch (option) {
        case 'p':
            graph1_fname = optarg;
            break;

        case 'q':
            graph2_fname = optarg;
            break;

        case 'd':
            delta = atof(optarg);
            break;

        default:
            cout << "invalid option: " << (char)optopt << endl;
            break;
        }
    }

    cout << "graph 1 filename: " << graph1_fname << endl;
    cout << "graph 2 filename: " << graph2_fname << endl;
    cout << "delta: " << delta << endl;

    return 0;
}
