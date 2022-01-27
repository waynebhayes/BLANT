#ifndef CDIJKSTRA_SEEDING_HEADER
#define CDIJKSTRA_SEEDING_HEADER

#include <string>
#include <utility>

#include "graph.h"

using std::pair;
using std::string;

namespace cdijkstra {

pair<unsigned int, unsigned int> get_seed_line(const string& seed_fname, unsigned int line_number, const Graph& g);

}

#endif
