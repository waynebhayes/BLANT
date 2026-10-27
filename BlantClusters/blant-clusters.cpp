#include "cluster.h"
#include "extsort.h"
#include "network.h"
#include "options.h"
#include "samples.h"

#include <algorithm>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

std::string g_tmpdir;

void removeTmpdir() {
    if (g_tmpdir.empty()) return;
    std::string cmd = "rm -rf '" + g_tmpdir + "'";
    if (std::system(cmd.c_str()) != 0) std::fprintf(stderr, "warning: could not remove %s\n",
                                                     g_tmpdir.c_str());
    g_tmpdir.clear();
}

void onSignal(int) {
    removeTmpdir();
    _exit(1);
}

std::string makeTmpdir(const std::string& parent) {
    std::string tmpl = parent + "/blant-clusters.XXXXXX";
    std::vector<char> buf(tmpl.begin(), tmpl.end());
    buf.push_back('\0');
    if (!mkdtemp(buf.data()))
        throw std::runtime_error("cannot create a temp dir under '" + parent + "'");
    return std::string(buf.data());
}

bool fileExists(const std::string& p) { return access(p.c_str(), R_OK) == 0; }

}  // namespace

int main(int argc, const char** argv) {
    Options opt;
    try {
        opt = parseOptions(argc, argv);
    } catch (const std::runtime_error& e) {
        std::fprintf(stderr, "FATAL ERROR: %s\n", e.what());
        return 1;
    }
    if (opt.ks.empty() && opt.netPath.empty()) return 1;  // usage was printed

    int status = 0;
    try {
        Network orig = loadNetwork(opt.netPath, opt.weighted);

        g_tmpdir = makeTmpdir(opt.tmpDir);
        std::atexit(removeTmpdir);
        std::signal(SIGINT, onSignal);
        std::signal(SIGTERM, onSignal);
        std::fprintf(stderr, "TMPDIR is %s\n", g_tmpdir.c_str());

        // BLANT is always handed the original network: -C makes it complement
        // internally. The clustering stages run on the complement.
        const std::string netForBlant = opt.netPath;
        Network net = orig;
        if (opt.complement) {
            const std::string cpath = g_tmpdir + (opt.weighted ? "/net.elw" : "/net.el");
            writeComplement(orig, cpath, opt.weighted);
            net = loadNetwork(cpath, opt.weighted);
        }

        // Run BLANT once per k, unless the caller supplied the output.
        std::vector<std::string> blantOut(opt.ks.size());
        for (size_t i = 0; i < opt.ks.size(); ++i) {
            const int k = opt.ks[i];
            if (!opt.blantDir.empty()) {
                blantOut[i] = opt.blantDir + "/blant" + std::to_string(k) + ".out";
                if (!fileExists(blantOut[i]))
                    throw std::runtime_error("no BLANT output for k=" + std::to_string(k) +
                                             " at '" + blantOut[i] + "'");
                continue;
            }
            const std::string maps = "canon_maps/canon_list" + std::to_string(k) + ".txt";
            if (!fileExists(maps)) {
                std::fprintf(stderr, "WARNING: no %s; skipping k=%d\n", maps.c_str(), k);
                blantOut[i].clear();
                continue;
            }
            blantOut[i] = g_tmpdir + "/blant" + std::to_string(k) + ".out";
            const int rc = runBlant(opt, k, blantOut[i], netForBlant);
            if (rc != 0) {
                std::fprintf(stderr, "FATAL ERROR: BLANT exited %d for k=%d\n", rc, k);
                status = 1;
                blantOut[i].clear();
            }
        }

        std::vector<Cluster> all;
        for (size_t e = 0; e < opt.eds.size(); ++e) {
            for (size_t i = 0; i < opt.ks.size(); ++i) {
                const int k = opt.ks[i];
                if (blantOut[i].empty()) continue;
                CanonMaps maps;
                if (!maps.load(k, "canon_maps")) {
                    std::fprintf(stderr, "WARNING: no canon_maps for k=%d; skipping\n", k);
                    continue;
                }
                ExternalSorter sorter(g_tmpdir + "/clus" + std::to_string(k) + "_" +
                                          std::to_string(e),
                                      opt.sortBufBytes);
                emitRecords(net, blantOut[i], maps, opt.communityMode, opt.eds[e], k,
                            opt.sqrNorm, opt.bestItem,
                            [&](const std::string& r) { sorter.add(r); });
                sorter.finish();

                Clusterer cl(net, opt, k, opt.eds[e]);
                std::string rec;
                long idx = 0;
                while (sorter.next(rec)) cl.addRecord(rec, ++idx);
                sorter.discard();

                std::vector<Cluster> clusters = cl.run();
                dedupExact(clusters);
                sortFileOrder(clusters);
                if (opt.onlyOne && !clusters.empty()) clusters.resize(1);
                for (auto& c : clusters) all.push_back(std::move(c));
            }
        }

        sortFileOrder(all);
        dedupOverlap(all, opt.overlap);
        std::stable_sort(all.begin(), all.end(), [](const Cluster& a, const Cluster& b) {
            if (a.nodes.size() != b.nodes.size()) return a.nodes.size() > b.nodes.size();
            return a.edges > b.edges;
        });
        for (auto& c : all) std::printf("%s\n", formatCluster(c, opt, net).c_str());
    } catch (const std::runtime_error& e) {
        std::fprintf(stderr, "FATAL ERROR: %s\n", e.what());
        return 1;
    }
    return status;
}
