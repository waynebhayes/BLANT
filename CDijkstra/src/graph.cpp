#include "graph.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include "matrix.h"

using std::ifstream;
using std::map;
using std::set;
using std::string;

using cdijkstra::Graph;
using cdijkstra::Matrix;

namespace {

unsigned int find_or_insert_node(map<string, unsigned int>& nodes, string node_name, unsigned int& next_node_index) {
    auto elem = nodes.find(node_name);

    if (elem != nodes.end()) {
        return elem->second;
    } else {
        unsigned int node = next_node_index++;
        nodes[node_name] = node;

        return node;
    }
}

}

Graph::Graph(string file)
 : adj_mat(0, 0, 0) {
    // open the file
    ifstream f(file);
    unsigned int n_nodes = 0;

    for (string line; getline(f, line);) {
        // split line by space
        unsigned int sep = line.find(' ');

        if (sep == string::npos) {
            // TODO: error, invalid format
        }

        // first node: 0 to sep
        // second node: sep + 1 to end
        string node1_name = line.substr(0, sep);
        string node2_name = line.substr(sep + 1);

        // get node indices, inserting into registry if required
        unsigned int node1 = find_or_insert_node(nodes, node1_name, n_nodes);
        unsigned int node2 = find_or_insert_node(nodes, node2_name, n_nodes);

        // add edge to edges
        auto node1_elem = edges.find(node1);
        if (node1_elem == edges.end()) {
            edges[node1] = set<unsigned int>();
        }

        edges[node1].insert(node2);
        
        auto node2_elem = edges.find(node2);
        if (node2_elem == edges.end()) {
            edges[node2] = set<unsigned int>();
        }

        edges[node2].insert(node1);
    }

    // alloc & init adjacency matrix
    adj_mat = Matrix<bool>(n_nodes, n_nodes, false);

    for (auto entry = edges.begin(); entry != edges.end(); ++entry) {
        unsigned int n = entry->first;

        for (auto m = entry->second.begin(); m != entry->second.end(); ++m) {
            adj_mat.at(n, *m) = true;
            adj_mat.at(*m, n) = true;
        }
    }
}

unsigned int Graph::size() const {
    return adj_mat.rows();
}

unsigned int Graph::index(const string& node) const {
    return nodes.at(node);
}
