#ifndef CDIJKSTRA_GRAPH_HEADER
#define CDIJKSTRA_GRAPH_HEADER

#include <map>
#include <set>
#include <string>

using std::map;
using std::set;
using std::string;

class Graph {
public:
    Graph(string file);
    ~Graph();

    Graph(const Graph& g) = delete;
    Graph& operator =(Graph& g) = delete;

    unsigned int size() const;
    unsigned int index(const string& node) const;
    // bool has_edge(int n, int m) const;
    // vector<int> neighbors(int n) const;

private:
    unsigned int n_nodes;
    bool** adjacency_matrix; // an adjacency matrix of size n_nodes by n_nodes
    map<unsigned int, set<unsigned int>> edges; // an adjacency list mapping each node to its edges
    map<string, unsigned int> nodes; // a registry of node name to int index in g
};

#endif
