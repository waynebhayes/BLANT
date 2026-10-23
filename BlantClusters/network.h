#pragma once
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

double choose2(long n);

// Throws std::runtime_error when the network cannot be clustered.
void checkpointNetworkStats(long n, long m);

struct Network {
    std::vector<std::string> names;                        // id -> name
    std::vector<std::vector<std::pair<int, double>>> adj;  // id -> (nbr, weight)
    bool weighted = false;
    long n = 0;
    long m = 0;
    double meanED = 0;
    std::unordered_map<std::string, int> index;

    // -1 when the name is not in the network.
    int idOf(const std::string& s) const;
};

// Reads an .el / .elw file. Blank lines and lines whose first non-space
// character is '#' are skipped. A self-loop is a fatal input error; a duplicate
// edge keeps the last weight and is counted once.
Network loadNetwork(const std::string& path, bool weighted);

// Streams every unordered non-edge of `net` to `path`. Two columns, or three
// with weight 1 when `weighted`, so the result loads back with the same flag.
void writeComplement(const Network& net, const std::string& path, bool weighted);
