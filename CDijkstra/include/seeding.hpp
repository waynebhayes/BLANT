#ifndef CDIJKSTRA_SEEDING_HEADER
#define CDIJKSTRA_SEEDING_HEADER

#include <string>
#include <utility>

#include "graph.h"
#include "matrix.h"

using std::pair;
using std::string;

using cdijkstra::Graph;
using cdijkstra::Matrix;

namespace cdijkstra {

pair<unsigned int, unsigned int> get_seed_line(const string& seed_fname, unsigned int line_number, const Graph& g);

}

#endif
