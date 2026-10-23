#include "network.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

double choose2(long n) { return n * (n - 1) / 2.0; }

void checkpointNetworkStats(long n, long m) {
    (void)m;
    if (n < 2)
        throw std::runtime_error("network has " + std::to_string(n) +
                                 " node(s); need at least 2 nodes to cluster");
}

int Network::idOf(const std::string& s) const {
    auto it = index.find(s);
    return it == index.end() ? -1 : it->second;
}

namespace {
unsigned long long edgeKey(int u, int v) {
    if (u > v) std::swap(u, v);
    return ((unsigned long long)(unsigned)u << 32) | (unsigned)v;
}
}  // namespace

Network loadNetwork(const std::string& path, bool weighted) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open network '" + path + "'");

    Network net;
    net.weighted = weighted;
    std::unordered_map<unsigned long long, std::pair<int, size_t>> seen;
    std::string line;
    long lineNo = 0;

    while (std::getline(in, line)) {
        ++lineNo;
        std::istringstream is(line);
        std::string a, b, w;
        if (!(is >> a)) continue;
        if (a[0] == '#') continue;
        if (!(is >> b)) break;                  // a lone token: nothing to add
        double weight = 1.0;
        if (is >> w) weight = std::strtod(w.c_str(), nullptr);
        if (a == b)
            throw std::runtime_error("self-loop at line " + std::to_string(lineNo) +
                                     " of '" + path + "'");

        auto add = [&](const std::string& name) {
            auto it = net.index.find(name);
            if (it != net.index.end()) return it->second;
            int id = (int)net.names.size();
            net.names.push_back(name);
            net.index[name] = id;
            net.adj.emplace_back();
            return id;
        };
        int u = add(a), v = add(b);

        unsigned long long key = edgeKey(u, v);
        auto it = seen.find(key);
        if (it != seen.end()) {
            // duplicate: last weight wins, counted once
            net.adj[it->second.first][it->second.second].second = weight;
            for (auto& e : net.adj[v])
                if (e.first == u) e.second = weight;
            continue;
        }
        seen[key] = {u, net.adj[u].size()};
        net.adj[u].emplace_back(v, weight);
        net.adj[v].emplace_back(u, weight);
        ++net.m;
    }

    net.n = (long)net.names.size();
    checkpointNetworkStats(net.n, net.m);
    net.meanED = net.m / choose2(net.n);
    return net;
}

void writeComplement(const Network& net, const std::string& path) {
    std::unordered_map<unsigned long long, char> present;
    for (int u = 0; u < net.n; ++u)
        for (auto& e : net.adj[u]) present[edgeKey(u, e.first)] = 1;

    std::string out;
    for (int u = 0; u < net.n; ++u)
        for (int v = u + 1; v < net.n; ++v)
            if (!present.count(edgeKey(u, v)))
                out += net.names[u] + "\t" + net.names[v] + "\t1\n";

    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write complement '" + path + "'");
    f << out;
}
