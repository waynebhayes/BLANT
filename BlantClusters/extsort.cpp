#include "extsort.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <stdexcept>
#include <vector>

double leadingNumber(const std::string& record) {
    return std::strtod(record.c_str(), nullptr);
}

namespace {
// Descending on the weight, then ascending on the whole record. The script's
// `sort -gr` breaks ties on the whole line descending, but its lines hold node
// names where these hold ids, so tie order cannot match it either way.
bool recordLess(const std::string& a, const std::string& b) {
    double wa = leadingNumber(a), wb = leadingNumber(b);
    if (wa != wb) return wa > wb;
    return a < b;
}
}  // namespace

struct ExternalSorter::Impl {
    std::string runPrefix;
    size_t bufBytes;
    std::vector<std::string> buf;
    size_t bufUsed = 0;
    long total = 0;
    std::vector<std::string> runs;

    struct Item {
        double weight;
        std::string record;
        size_t run;
        bool operator<(const Item& o) const {
            if (weight != o.weight) return weight < o.weight;  // max-heap
            return record > o.record;                          // ties: smaller first
        }
    };
    std::vector<std::unique_ptr<std::ifstream>> files;
    std::priority_queue<Item> heap;

    Impl(std::string prefix, size_t bytes) : runPrefix(std::move(prefix)), bufBytes(bytes) {}

    void spill() {
        if (buf.empty()) return;
        std::sort(buf.begin(), buf.end(), recordLess);
        std::string name = runPrefix + ".run" + std::to_string(runs.size());
        std::ofstream f(name, std::ios::binary);
        if (!f) throw std::runtime_error("cannot write sort run '" + name + "'");
        for (auto& r : buf) f << r << '\n';
        f.close();
        runs.push_back(name);
        buf.clear();
        bufUsed = 0;
    }

    void pushFrom(size_t run) {
        std::string line;
        if (std::getline(*files[run], line))
            heap.push({leadingNumber(line), std::move(line), run});
    }
};

ExternalSorter::ExternalSorter(std::string runPrefix, size_t bufBytes)
    : p(new Impl(std::move(runPrefix), bufBytes ? bufBytes : 1)) {}

ExternalSorter::~ExternalSorter() { discard(); }

void ExternalSorter::add(std::string record) {
    ++p->total;
    const size_t cost = sizeof(std::string) + record.size() + 1;
    if (!p->buf.empty() && p->bufUsed + cost > p->bufBytes) p->spill();
    p->bufUsed += cost;
    p->buf.push_back(std::move(record));
}

long ExternalSorter::size() const { return p->total; }

// shortcut: all runs are opened at once, so -B must leave the run count
// (total bytes / -B) under the open-file limit; add multi-pass merging if -B
// ever needs to be tiny.
void ExternalSorter::finish() {
    p->spill();
    p->files.clear();
    while (!p->heap.empty()) p->heap.pop();
    for (size_t r = 0; r < p->runs.size(); ++r) {
        p->files.emplace_back(new std::ifstream(p->runs[r], std::ios::binary));
        if (!*p->files.back())
            throw std::runtime_error("cannot read sort run '" + p->runs[r] + "'");
        p->pushFrom(r);
    }
}

bool ExternalSorter::next(std::string& out) {
    if (p->heap.empty()) return false;
    Impl::Item top = p->heap.top();
    p->heap.pop();
    out = std::move(top.record);
    p->pushFrom(top.run);
    return true;
}

void ExternalSorter::discard() {
    p->heap = std::priority_queue<Impl::Item>();
    p->files.clear();
    for (auto& r : p->runs) std::remove(r.c_str());
    p->runs.clear();
}
