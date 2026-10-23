// Independent invariant checker for blant-clusters output. It knows nothing
// about how the clusters were found; it only re-derives the facts the output
// claims, from the network.
//
// usage: check <network.el> <clusters.txt> <edgeDensity> <minEdgeMean> <overlap>
#include "../network.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

static unsigned long long key(int u, int v) {
    if (u > v) std::swap(u, v);
    return ((unsigned long long)(unsigned)u << 32) | (unsigned)v;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::fprintf(stderr, "usage: %s <network.el> <clusters.txt> <ED> <minEdgeMean> <overlap>\n",
                     argv[0]);
        return 2;
    }
    const double wantED = std::strtod(argv[3], nullptr);
    const double wantEdgeMean = std::strtod(argv[4], nullptr);
    const double overlap = std::strtod(argv[5], nullptr);

    Network net = loadNetwork(argv[1], false);
    std::unordered_set<unsigned long long> edges;
    for (int u = 0; u < net.n; ++u)
        for (auto& e : net.adj[(size_t)u]) edges.insert(key(u, e.first));

    FILE* in = std::fopen(argv[2], "r");
    if (!in) { std::perror(argv[2]); return 2; }

    std::vector<std::vector<int>> clusters;
    std::vector<long> claimedEdges;
    int bad = 0, lineNo = 0;
    char line[65536];
    while (std::fgets(line, sizeof line, in)) {
        ++lineNo;
        long nnodes = 0, ne = 0, maxE = 0;
        double dens = 0;
        int k = 0;
        if (std::sscanf(line, "%ld nodes, %ld/%ld edges, %lf%% density from k%d",
                        &nnodes, &ne, &maxE, &dens, &k) != 5) {
            std::fprintf(stderr, "line %d: cannot parse: %s", lineNo, line);
            ++bad;
            continue;
        }
        if (maxE != (long)(nnodes * (nnodes - 1) / 2)) {
            std::fprintf(stderr, "line %d: maxEdges %ld != C(%ld,2)\n", lineNo, maxE, nnodes);
            ++bad;
        }
        if (ne <= 0) { std::fprintf(stderr, "line %d: no edges\n", lineNo); ++bad; }
        if (maxE > 0 && (double)ne / (double)maxE < wantED - 1e-9) {
            std::fprintf(stderr, "line %d: density %g < %g\n", lineNo, (double)ne / maxE, wantED);
            ++bad;
        }

        const char* p = std::strstr(line, "nodeSet {");
        if (!p) { std::fprintf(stderr, "line %d: no nodeSet\n", lineNo); ++bad; continue; }
        p += 9;
        std::vector<int> nodes;
        std::set<int> uniq;
        std::string word;
        for (; *p && *p != '}'; ++p) {
            if (std::isspace((unsigned char)*p)) {
                if (!word.empty()) {
                    int id = net.idOf(word);
                    if (id < 0) {
                        std::fprintf(stderr, "line %d: node '%s' not in network\n",
                                     lineNo, word.c_str());
                        ++bad;
                    } else {
                        nodes.push_back(id);
                        uniq.insert(id);
                    }
                    word.clear();
                }
            } else {
                word += *p;
            }
        }
        if (uniq.size() != nodes.size()) {
            std::fprintf(stderr, "line %d: duplicate nodes in set\n", lineNo);
            ++bad;
        }
        if ((long)nodes.size() != nnodes) {
            std::fprintf(stderr, "line %d: %zu node names but claimed %ld\n",
                         lineNo, nodes.size(), nnodes);
            ++bad;
        }
        long induced = 0;
        double wsum = 0;
        for (size_t i = 0; i < nodes.size(); ++i)
            for (size_t j = i + 1; j < nodes.size(); ++j) {
                unsigned long long kk = key(nodes[i], nodes[j]);
                if (!edges.count(kk)) continue;
                ++induced;
                for (auto& e : net.adj[(size_t)nodes[i]])
                    if (e.first == nodes[j]) { wsum += e.second; break; }
            }
        if (induced != ne) {
            std::fprintf(stderr, "line %d: induced edges %ld but claimed %ld\n",
                         lineNo, induced, ne);
            ++bad;
        }
        if (ne > 0 && wsum / (double)ne < wantEdgeMean - 1e-9) {
            std::fprintf(stderr, "line %d: edgeMean %g < %g\n", lineNo, wsum / ne, wantEdgeMean);
            ++bad;
        }
        clusters.push_back(std::move(nodes));
        claimedEdges.push_back(ne);
    }
    std::fclose(in);

    if (clusters.empty()) {
        std::fprintf(stderr, "no clusters on stdout\n");
        return 1;
    }

    for (size_t i = 0; i < clusters.size(); ++i)
        for (size_t j = i + 1; j < clusters.size(); ++j) {
            size_t shared = 0;
            for (int u : clusters[i])
                if (std::find(clusters[j].begin(), clusters[j].end(), u) != clusters[j].end())
                    ++shared;
            const size_t smaller = std::min(clusters[i].size(), clusters[j].size());
            if ((double)shared > (double)smaller * overlap + 1e-9) {
                std::fprintf(stderr, "clusters %zu and %zu share %zu of %zu > %g\n",
                             i + 1, j + 1, shared, smaller, overlap);
                ++bad;
            }
        }

    if (bad) { std::fprintf(stderr, "%d invariant violation(s)\n", bad); return 1; }
    std::printf("check: %zu clusters ok\n", clusters.size());
    return 0;
}
