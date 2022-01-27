#ifndef CDIJKSTRA_SIM_HEADER
#define CDIJKSTRA_SIM_HEADER

#include <string>

#include "graph.h"
#include "matrix.h"

using std::string;

using cdijkstra::Graph;

namespace cdijkstra {

cdijkstra::Matrix<double> get_sim(const string& sim_file, const Graph& graph1, const Graph& graph2);

}

#endif
