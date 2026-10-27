#include "cluster.h"
#include "extsort.h"
#include "network.h"
#include "options.h"
#include "samples.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace {

// A plain buffer, not a std::string, because the signal handler reads it.
char g_tmpdir[PATH_MAX];

// fork/execlp/waitpid rather than system(): this also runs from the signal
// handler, and those calls are async-signal-safe where system() is not.
bool removeTree(const char* path) {
    const pid_t pid = fork();
    if (pid < 0) return false;
    if (pid == 0) {
        execlp("rm", "rm", "-rf", path, (char*)nullptr);
        _exit(127);
    }
    int st;
    while (waitpid(pid, &st, 0) < 0)
        if (errno != EINTR) return false;
    return WIFEXITED(st) && WEXITSTATUS(st) == 0;
}

void removeTmpdir() {
    if (!g_tmpdir[0]) return;
    if (!removeTree(g_tmpdir)) std::fprintf(stderr, "warning: could not remove %s\n", g_tmpdir);
    g_tmpdir[0] = '\0';
}

void onSignal(int) {
    if (g_tmpdir[0]) removeTree(g_tmpdir);
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
        if (opt.blantDir.empty() && access(opt.blantExe.c_str(), X_OK) != 0)
            throw std::runtime_error("BLANT executable '" + opt.blantExe +
                                     "' does not exist or is not executable");

        const std::string tmpdir = makeTmpdir(opt.tmpDir);
        if (tmpdir.size() >= sizeof g_tmpdir)
            throw std::runtime_error("temp dir path '" + tmpdir + "' is too long");
        std::strcpy(g_tmpdir, tmpdir.c_str());
        std::atexit(removeTmpdir);
        std::signal(SIGINT, onSignal);
        std::signal(SIGTERM, onSignal);
        std::fprintf(stderr, "TMPDIR is %s\n", tmpdir.c_str());

        // BLANT is always handed the original network: -C makes it complement
        // internally. The clustering stages run on the complement.
        const std::string netForBlant = opt.netPath;
        Network net = orig;
        if (opt.complement) {
            const std::string cpath = tmpdir + (opt.weighted ? "/net.elw" : "/net.el");
            writeComplement(orig, cpath, opt.weighted);
            net = loadNetwork(cpath, opt.weighted);
        }

        // Run BLANT once per k, unless the caller supplied the output. Like the
        // script, every k runs at once in its own child process (so peak memory
        // is the sum of the BLANT runs), and all are waited for before stage 1.
        std::vector<std::string> blantOut(opt.ks.size());
        std::vector<pid_t> blantPid(opt.ks.size(), -1);
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
            blantOut[i] = tmpdir + "/blant" + std::to_string(k) + ".out";
            std::fflush(nullptr);
            const pid_t pid = fork();
            if (pid < 0) throw std::runtime_error("cannot fork to run BLANT");
            if (pid == 0) {
                // The parent owns the temp dir and its cleanup.
                std::signal(SIGINT, SIG_DFL);
                std::signal(SIGTERM, SIG_DFL);
                int rc;
                try {
                    rc = runBlant(opt, k, blantOut[i], netForBlant);
                } catch (const std::runtime_error& e) {
                    std::fprintf(stderr, "FATAL ERROR: %s\n", e.what());
                    rc = 1;
                }
                std::fflush(nullptr);
                _exit(rc);
            }
            blantPid[i] = pid;
        }
        for (size_t i = 0; i < opt.ks.size(); ++i) {
            if (blantPid[i] < 0) continue;
            int st;
            while (waitpid(blantPid[i], &st, 0) < 0)
                if (errno != EINTR) throw std::runtime_error("lost track of a BLANT run");
            const int rc = WIFEXITED(st) ? WEXITSTATUS(st) : 128 + WTERMSIG(st);
            if (rc != 0) {
                std::fprintf(stderr, "FATAL ERROR: BLANT exited %d for k=%d\n", rc, opt.ks[i]);
                status = 1;
                blantOut[i].clear();
            }
        }

        const WeightIndex wgt = buildWeightIndex(net);

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
                ExternalSorter sorter(tmpdir + "/clus" + std::to_string(k) + "_" +
                                          std::to_string(e),
                                      opt.sortBufBytes);
                emitRecords(net, blantOut[i], maps, opt.communityMode, opt.eds[e], k,
                            opt.sqrNorm, opt.bestItem,
                            [&](const std::string& r) { sorter.add(r); });
                sorter.finish();

                Clusterer cl(net, wgt, opt, k, opt.eds[e]);
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
