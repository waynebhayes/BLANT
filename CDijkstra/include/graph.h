#ifndef CDIJKSTRA_GRAPH_HEADER
#define CDIJKSTRA_GRAPH_HEADER

#include <map>
#include <set>
#include <string>

#include "matrix.h"

using std::map;
using std::set;
using std::string;

using cdijkstra::Matrix;

namespace cdijkstra {

class Graph {
public:
    Graph(string file);

    unsigned int size() const;
    unsigned int index(const string& node) const;
    // bool has_edge(int n, int m) const;
    // vector<int> neighbors(int n) const;

private:
    Matrix<bool> adj_mat;
    map<unsigned int, set<unsigned int>> edges; // an adjacency list mapping each node to its edges
    map<string, unsigned int> nodes; // a registry of node name to int index in g
};

}

#endif
