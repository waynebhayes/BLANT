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
    std::vector<int> Slist;    // same members as S, in entry order
    std::deque<int> queue;
    // adjacency with weights, for O(1) "is there an edge" and weight lookups
    std::vector<std::unordered_map<int, double>> wgt;
    std::mt19937 rng;
};

// In-memory order used before dedupExact: size desc, edges desc, edgeSum desc,
// then ascending node ids.
void sortSubfinal(std::vector<Cluster>& clusters);

// File order for subfinal*.out: size desc, then first node id ascending.
void sortFileOrder(std::vector<Cluster>& clusters);

// How many of the first nodes entering S become "started". The script marks
// i = 1 .. while(i < size*(1-overlap)), i.e. ceil(size*(1-overlap)) - 1.
long startMarkCount(long clusterSize, double overlap);

// Sorts by sortSubfinal, then drops a cluster when an already-kept one is a
// subset of it. The sort is what guarantees a bigger set is never the one
// dropped.
void dedupExact(std::vector<Cluster>& clusters);

// Drops a cluster when it overlaps more than `overlap` of its own size with any
// already-kept one.
void dedupOverlap(std::vector<Cluster>& clusters, double overlap);

// One output line, no trailing newline.
std::string formatCluster(const Cluster& c, const Options& opt, const Network& net);
