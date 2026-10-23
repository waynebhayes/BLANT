#include "samples.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <sys/wait.h>
#include <unistd.h>

namespace {
std::string join(const std::vector<std::string>& v) {
    std::string out;
    for (auto& s : v) {
        if (s.empty()) continue;
        if (!out.empty()) out += ' ';
        out += s;
    }
    return out;
}
}  // namespace

bool CanonMaps::load(int k, const std::string& dir) {
    canonEdges.clear();
    orbit2canon.clear();
    const std::string cpath = dir + "/canon_list" + std::to_string(k) + ".txt";
    const std::string opath = dir + "/orbit_map" + std::to_string(k) + ".txt";
    std::ifstream c(cpath), o(opath);
    if (!c || !o) return false;

    std::string line;
    std::getline(c, line);  // count header
    while (std::getline(c, line)) {
        std::istringstream is(line);
        std::string f1, f2, f3;
        if (!(is >> f1 >> f2 >> f3)) { canonEdges.push_back(0); continue; }
        // field 3 of the raw file is the canon's edge count
        canonEdges.push_back(std::atoi(f3.c_str()));
    }

    std::getline(o, line);  // count header
    int maxOrbit = -1;
    std::vector<std::pair<int, int>> pending;
    for (long idx = 0; std::getline(o, line); ++idx) {
        std::istringstream is(line);
        for (int orbit; is >> orbit; ) {
            pending.emplace_back(orbit, (int)idx);
            maxOrbit = std::max(maxOrbit, orbit);
        }
    }
    orbit2canon.assign((size_t)maxOrbit + 1, 0);
    for (auto& p : pending) orbit2canon[(size_t)p.first] = p.second;
    return true;
}

int minEdgesFor(double edgeDensity, int k) {
    // The (int) cast reproduces the script's rounded_edC adjustment; the -1 is
    // unconditional there, and k-1 is the floor because a k-graphlet is connected.
    int rounded = (int)(edgeDensity * choose2(k));
    return std::max(rounded - 1, k - 1);
}

long emitRecords(const Network& net, const std::string& blantOut,
                 const CanonMaps& maps, char mode, double edgeDensity, int k,
                 bool sqrNorm, bool bestItem,
                 const std::function<void(const std::string&)>& emit) {
    std::ifstream in(blantOut);
    if (!in) throw std::runtime_error("cannot open BLANT output '" + blantOut + "'");

    const int minEdges = minEdgesFor(edgeDensity, k);
    std::map<int, long> kc;            // item -> count, ordered so output is stable
    std::map<int, std::vector<int>> nbrs;
    long total = 0, records = 0, lineNo = 0;
    int cur = -1;
    bool first = true;
    std::string line;

    auto flush = [&]() {
        if (kc.empty()) return;
        std::vector<int> items;
        if (bestItem) {
            int best = kc.begin()->first;
            for (auto& kv : kc)
                if (kv.second > kc[best]) best = kv.first;
            items.push_back(best);
        } else {
            for (auto& kv : kc) items.push_back(kv.first);
        }
        for (int item : items) {
            char head[128];
            double w = sqrNorm ? (double)kc[item] * (double)kc[item] / (double)total
                               : (double)kc[item];
            std::snprintf(head, sizeof head, "%.17g %d %d", w, cur, item);
            std::string rec = head;
            for (int v : nbrs[item]) rec += " " + std::to_string(v);
            emit(rec);
            ++records;
        }
        kc.clear();
        nbrs.clear();
        total = 0;
    };

    while (std::getline(in, line)) {
        ++lineNo;
        std::istringstream is(line);
        std::string name;
        int itemRaw = 0;
        long count = 0;
        if (!(is >> name >> itemRaw >> count)) continue;
        if (count <= 0) continue;

        int node = net.idOf(name);
        if (node < 0)
            throw std::runtime_error("BLANT output line " + std::to_string(lineNo) +
                                     " names node '" + name + "' which is not in the network");
        if (first) { cur = node; first = false; }
        else if (node < cur)
            throw std::runtime_error("BLANT output nodes are not ascending at line " +
                                     std::to_string(lineNo) + " ('" + name + "')");
        else if (node > cur) { flush(); cur = node; }

        int item = itemRaw;
        int edges;
        if (mode == 'o') {
            if (itemRaw < 0 || (size_t)itemRaw >= maps.orbit2canon.size())
                throw std::runtime_error("orbit " + std::to_string(itemRaw) +
                                         " has no canon in orbit_map" + std::to_string(k));
            edges = maps.canonEdges[(size_t)maps.orbit2canon[(size_t)itemRaw]];
        } else {
            if (itemRaw < 0 || (size_t)itemRaw >= maps.canonEdges.size())
                throw std::runtime_error("canon " + std::to_string(itemRaw) +
                                         " is not in canon_list" + std::to_string(k));
            edges = maps.canonEdges[(size_t)itemRaw];
        }
        if (edges < minEdges) continue;

        kc[item] += count;
        total += count;
        for (std::string v; is >> v; ) {
            int id = net.idOf(v);
            if (id < 0)
                throw std::runtime_error("BLANT output line " + std::to_string(lineNo) +
                                         " names neighbour '" + v + "' which is not in the network");
            nbrs[item].push_back(id);
        }
    }
    flush();
    return records;
}

std::string blantCommand(const Options& opt, int k, const std::string& netForBlant,
                         const std::string& turanCount) {
    std::vector<std::string> tok;
    tok.push_back(opt.blantCmd);
    tok.push_back(opt.precision);
    if (opt.quiet) tok.push_back("-q");
    if (opt.weighted) tok.push_back("-w");
    if (opt.turan && !turanCount.empty()) tok.push_back("-A " + turanCount);
    if (opt.complement) tok.push_back("-C");
    tok.push_back(opt.seedArg);
    tok.push_back("-k" + std::to_string(k));
    tok.push_back(opt.sampleMethod);
    tok.push_back(std::string("-mc") + opt.communityMode);
    tok.push_back(netForBlant);
    return join(tok);
}

int runBlant(const Options& opt, int k, const std::string& outPath,
             const std::string& netForBlant) {
    if (access(opt.blantExe.c_str(), X_OK) != 0)
        throw std::runtime_error("BLANT executable '" + opt.blantExe +
                                 "' does not exist or is not executable");
    std::string turanCount;
    if (opt.turan) {
        std::string probe = "./scripts/absolute-clique-count.sh " + std::to_string(k) +
                            " '" + netForBlant + "' 2>/dev/null";
        FILE* pipe = popen(probe.c_str(), "r");
        if (!pipe)
            throw std::runtime_error("cannot run scripts/absolute-clique-count.sh for -A");
        char buf[128] = {0};
        if (fgets(buf, sizeof buf, pipe)) turanCount = buf;
        while (!turanCount.empty() && isspace((unsigned char)turanCount.back()))
            turanCount.pop_back();
        pclose(pipe);
        if (turanCount.empty())
            throw std::runtime_error("-A needs ./scripts/absolute-clique-count.sh, which did not"
                                     " produce a count for k=" + std::to_string(k));
    }
    std::string cmd = blantCommand(opt, k, netForBlant, turanCount) + " > '" + outPath + "'";
    int rc = std::system(cmd.c_str());
    if (rc == -1) throw std::runtime_error("failed to run '" + cmd + "'");
    return WEXITSTATUS(rc);
}
