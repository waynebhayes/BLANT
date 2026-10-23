#include "cluster.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

Clusterer::Clusterer(const Network& n, const Options& o, int kk, double ed)
    : net(n), opt(o), k(kk), edgeDensity(ed),
      count((size_t)n.n, 0.0), hasCount((size_t)n.n, 0), lineOf((size_t)n.n, -1),
      gnbr((size_t)n.n), inS((size_t)n.n, 0), visitedQ((size_t)n.n, 0),
      rng(o.seed ? (unsigned)o.seed : std::random_device{}()) {}

void Clusterer::addRecord(const std::string& record, long recordIndex) {
    std::istringstream is(record);
    double w;
    int u, item;
    if (!(is >> w >> u >> item))
        throw std::runtime_error("malformed stage-1 record '" + record + "'");
    if (u < 0 || u >= net.n)
        throw std::runtime_error("stage-1 record names node " + std::to_string(u) +
                                 " outside a network of " + std::to_string(net.n));
    count[(size_t)u] += w;
    hasCount[(size_t)u] = 1;
    lineOf[(size_t)u] = recordIndex;
    if ((long)nodeAt.size() <= recordIndex) nodeAt.resize(recordIndex + 1, -1);
    nodeAt[(size_t)recordIndex] = u;
    recordCount = recordIndex;

    const double thresh = opt.minEdgeMean / 2;
    for (int v; is >> v; ) {
        if (v < 0 || v >= net.n) continue;
        const bool linked = std::any_of(net.adj[(size_t)u].begin(), net.adj[(size_t)u].end(),
                                        [&](const std::pair<int, double>& e) {
                                            return e.first == v && e.second > thresh;
                                        });
        if (linked) {
            gnbr[(size_t)u][v] = 1;
            gnbr[(size_t)v][u] = 1;
        }
    }
}

bool Clusterer::highRelClusCount(int u, int v) const {
    if (!hasCount[(size_t)u] || count[(size_t)u] == 0) return true;
    if (!hasCount[(size_t)v] || count[(size_t)v] == 0) return false;
    return count[(size_t)v] / count[(size_t)u] >= std::sqrt(edgeDensity) * 0.9;
}

int Clusterer::edgesIntoS(int v) const {
    int hits = 0;
    for (int u = 0; u < net.n; ++u) {
        if (!inS[(size_t)u] || u == v) continue;
        for (auto& e : net.adj[(size_t)u])
            if (e.first == v) { ++hits; break; }
    }
    return hits;
}

double Clusterer::wgtEdgesIntoS(int v) const {
    double sum = 0;
    for (int u = 0; u < net.n; ++u) {
        if (!inS[(size_t)u] || u == v) continue;
        for (auto& e : net.adj[(size_t)u])
            if (e.first == v) { sum += e.second; break; }
    }
    return sum;
}

void Clusterer::inDegrees(std::vector<std::pair<int, int>>& out) const {
    out.clear();
    for (int u = 0; u < net.n; ++u) {
        if (!inS[(size_t)u]) continue;
        int d = 0;
        for (auto& e : net.adj[(size_t)u])
            if (inS[(size_t)e.first]) ++d;
        out.emplace_back(u, d);
    }
}

void Clusterer::appendNeighbors(int u, int origin) {
    std::vector<int> nbrs;
    for (auto& kv : gnbr[(size_t)u]) nbrs.push_back(kv.first);
    std::shuffle(nbrs.begin(), nbrs.end(), rng);   // the script's randsort
    const long originLine = lineOf[(size_t)origin];
    for (int v : nbrs) {
        if (visitedQ[(size_t)v]) continue;
        if (lineOf[(size_t)v] < 0 || lineOf[(size_t)v] <= originLine) continue;
        if (!highRelClusCount(u, v)) continue;
        visitedQ[(size_t)v] = 1;
        queue.push_back(v);
    }
}

std::vector<Cluster> Clusterer::run() {
    long minClus;
    if (opt.minClusArg >= 1) {
        minClus = (long)opt.minClusArg;
    } else {
        const double eps = net.m / choose2(net.n);
        const double pVal = opt.minClusArg > 0 ? opt.minClusArg
                                              : 1.0 / (choose2(net.n) * choose2(net.n));
        long e = 1;
        for (; e < net.m; ++e)
            if (std::pow(eps, (double)e) < pVal) break;
        for (minClus = 2; minClus < net.n; ++minClus)
            if (edgeDensity * choose2(minClus) >= (double)e) break;
        if (opt.verbose)
            std::fprintf(stderr, "minClus %ld (%ld edges = pVal %g < %g for ED %g)\n",
                         minClus, e, std::pow(eps, (double)e), pVal, edgeDensity);
    }
    if (minClus < 1) minClus = 1;

    std::vector<Cluster> accepted;
    std::vector<char> started((size_t)net.n, 0);
    std::vector<int> Slist;
    std::vector<std::pair<int, int>> deg;
    const long startCount = (long)((double)recordCount * edgeDensity);

    for (long start = 1; start <= startCount; ++start) {
        if ((size_t)start >= nodeAt.size()) break;
        const int origin = nodeAt[(size_t)start];
        if (origin < 0 || started[(size_t)origin]) continue;
        started[(size_t)origin] = 1;

        std::fill(inS.begin(), inS.end(), 0);
        std::fill(visitedQ.begin(), visitedQ.end(), 0);
        Slist.clear();
        Sorder.clear();
        queue.clear();
        long misses = 0;
        long edgeCount = 0;
        double wgtEdgeCount = 0;

        visitedQ[(size_t)origin] = 1;
        queue.push_back(origin);

        while (!queue.empty()) {
            const int u = queue.front();
            queue.pop_front();
            const int newEdgeHits = edgesIntoS(u);
            const double wgtHits = wgtEdgesIntoS(u);
            edgeCount += newEdgeHits;
            wgtEdgeCount += wgtHits;
            inS[(size_t)u] = 1;
            Slist.push_back(u);
            const long Slen = (long)Slist.size();
            if ((double)Slen >= opt.maxClus) break;
            Sorder.push_back(u);

            if (Slen == 1) {
                appendNeighbors(u, origin);
            } else if ((double)edgeCount / choose2(Slen) < edgeDensity * 0.95) {
                inS[(size_t)u] = 0;              // u drops the density too low
                Slist.pop_back();
                Sorder.pop_back();
                edgeCount -= newEdgeHits;
                wgtEdgeCount -= wgtHits;
                if (++misses > std::max(Slen, net.n / 100)) break;
            } else {
                misses = 0;
                appendNeighbors(u, origin);
            }
        }

        // post-process: drop members with a very weak hold on the cluster
        if (queue.empty() && Slist.size() > 1) {
            inDegrees(deg);
            double sum = 0, sum2 = 0;
            std::unordered_map<int, int> freq;
            for (auto& p : deg) { sum += p.second; sum2 += (double)p.second * p.second;
                                  ++freq[p.second]; }
            const double nd = (double)deg.size();
            const double mean = sum / nd;
            const double sd = nd > 1 ? std::sqrt(std::fabs(sum2 - sum * sum / nd) / (nd - 1)) : 0;
            int maxFreq = 0, degMode = 0;
            for (auto& kv : freq)
                if (kv.second > maxFreq || (kv.second == maxFreq && kv.first < degMode)) {
                    maxFreq = kv.second; degMode = kv.first;
                }
            std::sort(deg.begin(), deg.end(),
                      [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                          return a.second != b.second ? a.second < b.second : a.first < b.first;
                      });
            for (auto& p : deg) {
                if (!inS[(size_t)p.first]) continue;
                if (p.second < mean - 3 * sd || p.second < degMode / 3.0) {
                    inS[(size_t)p.first] = 0;
                    edgeCount -= edgesIntoS(p.first);
                    wgtEdgeCount -= wgtEdgesIntoS(p.first);
                    auto it = std::find(Slist.begin(), Slist.end(), p.first);
                    if (it != Slist.end()) Slist.erase(it);
                }
            }
        }

        while (Slist.size() >= 2 &&
               (double)edgeCount / choose2((long)Slist.size()) < edgeDensity) {
            inDegrees(deg);
            int worst = -1, worstDeg = 0;
            for (auto& p : deg)
                if (worst < 0 || p.second < worstDeg) { worst = p.first; worstDeg = p.second; }
            if (worst < 0) break;
            inS[(size_t)worst] = 0;
            edgeCount -= edgesIntoS(worst);
            wgtEdgeCount -= wgtEdgesIntoS(worst);
            auto it = std::find(Slist.begin(), Slist.end(), worst);
            if (it != Slist.end()) Slist.erase(it);
        }

        const long Slen = (long)Slist.size();
        if (Slen >= minClus &&
            (edgeCount == 0 ? 0.0 : wgtEdgeCount / (double)edgeCount) >= opt.minEdgeMean) {
            Cluster c;
            c.nodes = Slist;
            std::sort(c.nodes.begin(), c.nodes.end());
            c.edges = edgeCount;
            c.edgeSum = wgtEdgeCount;
            c.k = k;
            accepted.push_back(std::move(c));
            const long mark = (long)((double)Slen * (1.0 - opt.overlap));
            for (long i = 0; i < mark && (size_t)i < Sorder.size(); ++i)
                started[(size_t)Sorder[(size_t)i]] = 1;
        }
    }
    return accepted;
}
