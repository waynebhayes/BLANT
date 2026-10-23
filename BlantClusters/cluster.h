#pragma once
#include "network.h"
#include "options.h"
#include <deque>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

struct Cluster {
    std::vector<int> nodes;   // ascending ids
    long edges = 0;
    double edgeSum = 0;
    int k = 0;
};

// Grows one dense community per start node, the way the script's second awk
// does: a greedy BFS over the graphlet-neighbour graph, a low-in-cluster-degree
// prune, then a peel until the edge density reaches the target.
class Clusterer {
public:
    Clusterer(const Network& net, const Options& opt, int k, double edgeDensity);

    // Stage-2 ingestion, called once per sorted record, in order.
    void addRecord(const std::string& record, long recordIndex);

    // Runs the start loop and returns the accepted clusters, duplicates and all.
    std::vector<Cluster> run();

private:
    void appendNeighbors(int u, int origin);
    bool highRelClusCount(int u, int v) const;
    int edgesIntoS(int v) const;
    double wgtEdgesIntoS(int v) const;
    // Fills `out` with (node, in-S degree) for everything currently in S.
    void inDegrees(std::vector<std::pair<int, int>>& out) const;

    const Network& net;
    const Options& opt;
    int k;
    double edgeDensity;
    long recordCount = 0;

    std::vector<double> count;
    std::vector<char> hasCount;
    std::vector<long> lineOf;
    std::vector<int> nodeAt;   // recordIndex -> node id (1-based)
    std::vector<std::unordered_map<int, char>> gnbr;
    std::vector<char> inS;
    std::vector<char> visitedQ;
    std::vector<int> Sorder;   // the order nodes entered S
    std::deque<int> queue;
    std::mt19937 rng;
};
