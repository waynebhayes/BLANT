#include "sim.h"

#include <fstream>
#include <iostream>
#include <string>

#include "graph.h"
#include "matrix.h"

using std::ifstream;
using std::string;

using cdijkstra::Graph;
using cdijkstra::Matrix;

namespace cdijkstra {

Matrix<double> get_sim(const string& sim_fname, const Graph& graph1, const Graph& graph2) {
    Matrix<double> sim(graph1.size(), graph2.size(), 0.0);

    // open sim file
    ifstream sim_file(sim_fname);

    // for each row:
    for (string line; std::getline(sim_file, line);) {
        unsigned int sep1 = line.find(' ');
        unsigned int sep2 = line.rfind(' ');

        if (sep1 == string::npos || sep2 == string::npos) {
            // TODO: invalid format, error
        }

        // first node: 0 to sep1
        // second node: sep1 + 1 to sep2
        // sim value: sep2 + 1 to end
        string g1_node = line.substr(0, sep1);
        string g2_node = line.substr(sep1 + 1, sep2 - sep1 - 1);
        double sim_val = stof(line.substr(sep2 + 1));

        int n = graph1.index(g1_node);
        int m = graph2.index(g2_node);
        sim.at(n, m) = sim_val;
    }

    return sim;
}

}
