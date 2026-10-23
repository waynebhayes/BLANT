#include "network.h"
#include "options.h"
#include <cstdio>

int main(int argc, const char** argv) {
    Options opt = parseOptions(argc, argv);
    if (opt.ks.empty() && opt.netPath.empty()) return 1;
    Network net = loadNetwork(opt.netPath, opt.weighted);
    std::fprintf(stderr, "%ld nodes, %ld edges, meanED %g\n", net.n, net.m, net.meanED);
    return 0;
}
