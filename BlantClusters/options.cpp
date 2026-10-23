#include "options.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <functional>
#include <sstream>
#include <stdexcept>

static const char* USAGE =
    "USAGE: blant-clusters [OPTIONS] blant_cmd Ks EDs network.el\n"
    "PURPOSE: use random samples of k-graphlets from BLANT in attempt to find large communities in network.el.\n"
    "REQUIRED ARGUMENTS:\n"
    "    blant_cmd is the raw BLANT command; usually it's just './blant', but can include options (eg './blant -C')\n"
    "    Ks is a list of graphlet sizes to do BLANT sampling in\n"
    "    EDs is a list of the edge density thresholds to explore\n"
    "OPTIONS (added BEFORE the blant.exe name)\n"
    "    -D DIR: use existing blant output files in directory DIR\n"
    "    -T DIR: use temp space here rather than /tmp\n"
    "    -A: estimate absolute counts with Turan's Shadow\n"
    "    -C: run on the complement graph, G'\n"
    "    -1: exit after printing only one cluster per k/ED\n"
    "    -h: use only the highest count graphlet/orbit for neighbors\n"
    "    -S: cluster count sorted by the normalized square\n"
    "    -w: networks are edge-weighted; pass -w to BLANT to use weighted graphlets\n"
    "    -v: verbose output\n"
    "    -nm: 'no members': suppress list of nodes in the community\n"
    "    -pD: tell BLANT to estimate graphlet counts with D digits of precision (default 1)\n"
    "    -nN: tell BLANT to take exactly N samples (mutually exclusive with -p above)\n"
    "    -r SEED: use the integer SEED as the random seed\n"
    "    -emt value: minimum edge mean threshold [default 0.5]\n"
    "    -m smallest: minimum cluster size, or a p-value when < 1\n"
    "    -M largest: stop building a community at this size (default: no limit)\n"
    "    -sSAMPLE_METHOD: BLANT's sampling method (MCMC, NBE, EBE, or RES)\n"
    "    -o OVERLAP: maximum fractional overlap with any one previous cluster (default 0.5)\n"
    "    -B BYTES: sort buffer budget in bytes (default 536870912)\n"
    "    -P N: accepted and ignored (edge prediction is not implemented)\n";

void printUsage(FILE* out) { std::fputs(USAGE, out); }

namespace {

std::string nextArg(int argc, const char* const* argv, int& i, const char* opt) {
    if (i + 1 >= argc)
        throw std::runtime_error(std::string("option '") + opt + "' needs an argument");
    return argv[++i];
}

double toDouble(const std::string& s, const char* opt) {
    char* end = nullptr;
    double d = std::strtod(s.c_str(), &end);
    if (end == s.c_str() || (end && *end))
        throw std::runtime_error(std::string("option '") + opt + "' got a non-number '" + s + "'");
    return d;
}

long toLong(const std::string& s, const char* opt) {
    char* end = nullptr;
    long v = std::strtol(s.c_str(), &end, 10);
    if (end == s.c_str() || (end && *end))
        throw std::runtime_error(std::string("option '") + opt + "' got a non-integer '" + s + "'");
    return v;
}

bool startsWith(const std::string& s, const char* p) {
    return s.compare(0, std::char_traits<char>::length(p), p) == 0;
}

}  // namespace

Options parseOptions(int argc, const char* const* argv) {
    Options o;
    if (argc == 1) {
        printUsage(stderr);
        return o;
    }

    int i = 1;
    // Longest prefixes first: -nm must be tested before -n.
    while (i < argc && argv[i][0] == '-') {
        std::string a = argv[i];
        if (a == "-1") o.onlyOne = true;
        else if (a == "-h") o.bestItem = true;
        else if (a == "-S") o.sqrNorm = true;
        else if (a == "-w") { o.weighted = true; }
        else if (a == "-A") o.turan = true;
        else if (a == "-C") o.complement = true;
        else if (a == "-v") { o.verbose = true; o.quiet = false; }
        else if (a == "-nm") o.printMembers = false;
        else if (a == "-emt") o.minEdgeMean = toDouble(nextArg(argc, argv, i, "-emt"), "-emt");
        else if (a == "-D") o.blantDir = nextArg(argc, argv, i, "-D");
        else if (a == "-T") o.tmpDir = nextArg(argc, argv, i, "-T");
        else if (a == "-B") o.sortBufBytes = (size_t)toLong(nextArg(argc, argv, i, "-B"), "-B");
        else if (a == "-s") o.sampleMethod = "-s" + nextArg(argc, argv, i, "-s");
        else if (startsWith(a, "-s")) o.sampleMethod = a;
        else if (a == "-n") o.precision = "-n " + nextArg(argc, argv, i, "-n");
        else if (startsWith(a, "-n") && a.size() > 2 && isdigit((unsigned char)a[2])) o.precision = a;
        else if (startsWith(a, "-p") && a.size() > 2 && isdigit((unsigned char)a[2])) o.precision = a + "L";
        else if (a == "-r") o.seedArg = "-r " + nextArg(argc, argv, i, "-r");
        else if (startsWith(a, "-r") && a.size() > 2) o.seedArg = a;
        else if (a == "-m") o.minClusArg = toDouble(nextArg(argc, argv, i, "-m"), "-m");
        else if (startsWith(a, "-m") && a.size() > 2) o.minClusArg = toDouble(a.substr(2), "-m");
        else if (a == "-M") o.maxClus = toDouble(nextArg(argc, argv, i, "-M"), "-M");
        else if (startsWith(a, "-M") && a.size() > 2) o.maxClus = toDouble(a.substr(2), "-M");
        else if (a == "-o") o.overlap = toDouble(nextArg(argc, argv, i, "-o"), "-o");
        else if (startsWith(a, "-o") && a.size() > 2) o.overlap = toDouble(a.substr(2), "-o");
        else if (startsWith(a, "-P")) { /* accepted, ignored */ }
        else throw std::runtime_error("unknown option '" + a + "'");
        ++i;
    }

    int rest = argc - i;
    if (rest != 4)
        throw std::runtime_error("expecting exactly 4 mandatory arguments, but you supplied " +
                                 std::to_string(rest));

    o.blantCmd = argv[i];
    o.blantExe = o.blantCmd.substr(0, o.blantCmd.find(' '));

    std::istringstream ks(argv[i + 1]);
    for (int k; ks >> k; ) o.ks.push_back(k);
    std::sort(o.ks.begin(), o.ks.end(), std::greater<int>());

    std::istringstream eds(argv[i + 2]);
    for (std::string t; eds >> t; ) o.eds.push_back(toDouble(t, "EDs"));

    o.netPath = argv[i + 3];

    if (o.seedArg.empty()) o.seed = 0;
    else {
        std::string s = o.seedArg;
        if (startsWith(s, "-r")) s = s.substr(2);
        o.seed = s.empty() ? 0 : toLong(s, "-r");
    }

    if (o.tmpDir.empty()) o.tmpDir = "/tmp";
    return o;
}
