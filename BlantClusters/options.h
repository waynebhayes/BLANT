#pragma once
#include <cstddef>
#include <cstdio>
#include <string>
#include <vector>

struct Options {
    std::string blantCmd;            // positional 1, e.g. "./blant -C"
    std::string blantExe;            // first whitespace-separated token of blantCmd
    std::vector<int> ks;             // positional 2, sorted descending
    std::vector<double> eds;         // positional 3
    std::string netPath;             // positional 4
    std::string blantDir;            // -D, empty means "run BLANT"
    std::string tmpDir;              // -T, default "/tmp"
    bool turan = false;              // -A
    bool complement = false;         // -C
    bool onlyOne = false;            // -1
    bool bestItem = false;           // -h
    bool sqrNorm = false;            // -S
    bool weighted = false;           // -w
    bool verbose = false;            // -v
    bool printMembers = true;        // -nm clears it
    bool quiet = true;               // -v clears it
    std::string precision = "-p0.9L";// -p0.9 / -p, -n
    std::string sampleMethod = "-sEBE"; // -s
    std::string seedArg;             // "-r N", empty when unset
    long seed = 0;                   // parsed from seedArg, 0 means entropy
    double minEdgeMean = 0.5;        // -emt
    double minClusArg = 0;           // -m
    double maxClus = 2147483647.0;   // -M
    double overlap = 0.5;            // -o
    char communityMode = 'g';        // the script hardcodes -mcg; no flag for 'o'
    size_t sortBufBytes = 512u << 20;// -B
};

void printUsage(FILE* out);

// Throws std::runtime_error on unknown options, a wrong positional argument
// count, or a bad numeric argument. With argc == 1 prints the usage text and
// returns a default-constructed Options (the caller decides what to do).
Options parseOptions(int argc, const char* const* argv);
