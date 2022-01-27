#include "seeding.h"

#include <fstream>
#include <string>
#include <utility>

using std::getline;
using std::ifstream;
using std::pair;
using std::string;

namespace cdijkstra {

pair<unsigned int, unsigned int> get_seed_line(const string& seed_fname, unsigned int line_number, const Graph& g) {
    // open file
    ifstream seed_file(seed_fname);
    string line;
    
    // skip to line
    for (unsigned int i = 0; i < line_number - 1; i++) {
        line.clear();
        getline(seed_file, line);
    }

    // read line, strip into (kval, node name)
    unsigned int sep = line.find(' ');

    if (sep == string::npos) {
        // TODO: error, invalid format
    }

    // return pair(kval, g.nodes[node])
    unsigned int kval = stoi(line.substr(0, sep));
    string node_name = line.substr(sep + 1);
    unsigned int node = g.index(node_name);

    return pair<unsigned int, unsigned int>(kval, node);
}

}
