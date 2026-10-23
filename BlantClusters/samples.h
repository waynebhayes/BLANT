#pragma once
#include "network.h"
#include "options.h"
#include <functional>
#include <string>
#include <vector>

struct CanonMaps {
    std::vector<int> canonEdges;    // canon -> edge count
    std::vector<int> orbit2canon;   // orbit -> canon

    // Reads <dir>/canon_list<k>.txt and <dir>/orbit_map<k>.txt.
    // Returns false (leaving the maps empty) if either file is missing.
    bool load(int k, const std::string& dir);
};

// Minimum edge count for a canon to be considered at this edge density.
int minEdgesFor(double edgeDensity, int k);

// Streams `blantOut` and hands one record per (node, item) to `emit`, as
// "<weight> <nodeId> <item> <neighbour>...". BLANT writes nodes in ascending
// order with the canons for a node contiguous, so this is a one-node-at-a-time
// accumulator. Throws std::runtime_error if the node ids are not ascending or
// a printed node is not in the network.
long emitRecords(const Network& net, const std::string& blantOut,
                 const CanonMaps& maps, char mode, double edgeDensity, int k,
                 bool sqrNorm, bool bestItem,
                 const std::function<void(const std::string&)>& emit);

// The BLANT command line, without the stdout redirection. `turanCount` is the
// absolute clique count to pass with -A; ignored unless opt.turan is set.
std::string blantCommand(const Options& opt, int k, const std::string& netForBlant,
                         const std::string& turanCount = "");

// Runs BLANT with stdout redirected to `outPath`; returns its exit status.
// Throws std::runtime_error if the executable is not runnable.
int runBlant(const Options& opt, int k, const std::string& outPath,
             const std::string& netForBlant);
