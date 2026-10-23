#include "tap.h"
#include "../options.h"
#include "../network.h"
#include "../extsort.h"
#include <cstring>
#include <stdexcept>

static Options parse(const char* const* argv, int argc) {
    return parseOptions(argc, argv);
}

static void testOptions() {
    { const char* a[] = {"./blant-clusters", "./blant", "3 5 4", "0.5 0.9", "x.el"};
      Options o = parse(a, 5);
      CHECK_EQ(o.blantExe, std::string("./blant"));
      CHECK_EQ(o.ks.size(), size_t(3));
      CHECK_EQ(o.ks[0], 5); CHECK_EQ(o.ks[1], 4); CHECK_EQ(o.ks[2], 3); // descending
      CHECK_EQ(o.eds.size(), size_t(2));
      CHECK_EQ(o.eds[1], 0.9);
      CHECK_EQ(o.netPath, std::string("x.el")); }

    { const char* a[] = {"./blant-clusters", "-m5", "-o", "0.7", "-emt", "0.6",
                         "-M", "12", "-r", "42", "-p0.9", "-w", "-v", "-nm",
                         "-1", "-h", "-S", "-A", "-C", "./blant -q", "4", "0.5", "x.elw"};
      Options o = parse(a, 23);
      CHECK_EQ(o.minClusArg, 5.0);
      CHECK_EQ(o.overlap, 0.7);
      CHECK_EQ(o.minEdgeMean, 0.6);
      CHECK_EQ(o.maxClus, 12.0);
      CHECK_EQ(o.seed, 42L);
      CHECK_EQ(o.precision, std::string("-p0.9L"));
      CHECK(o.weighted); CHECK(o.verbose); CHECK(!o.printMembers);
      CHECK(o.onlyOne); CHECK(o.bestItem); CHECK(o.sqrNorm);
      CHECK(o.turan); CHECK(o.complement);
      CHECK_EQ(o.blantExe, std::string("./blant"));
      CHECK_EQ(o.blantCmd, std::string("./blant -q")); }

    { const char* a[] = {"./blant-clusters", "-T", "/scratch", "-D", "out",
                         "-B", "1024", "./blant", "5", "0.5", "x.el"};
      Options o = parse(a, 11);
      CHECK_EQ(o.tmpDir, std::string("/scratch"));
      CHECK_EQ(o.blantDir, std::string("out"));
      CHECK_EQ(o.sortBufBytes, size_t(1024)); }

    // wrong arity
    { const char* a[] = {"./blant-clusters", "./blant", "5", "x.el"};
      bool threw = false;
      try { parse(a, 4); } catch (const std::runtime_error&) { threw = true; }
      CHECK(threw); }

    // unknown option
    { const char* a[] = {"./blant-clusters", "-Z", "./blant", "5", "0.5", "x.el"};
      bool threw = false;
      try { parse(a, 6); } catch (const std::runtime_error& e) {
          threw = std::strstr(e.what(), "unknown option") != nullptr; }
      CHECK(threw); }
}

static void testNetwork() {
    { Network n = loadNetwork("test/fixture/fixture.el", false);
      CHECK_EQ(n.n, 8L);
      CHECK_EQ(n.m, 13L);
      CHECK_EQ(n.names[0], std::string("a"));
      CHECK_EQ(n.idOf("h"), 7);
      CHECK_EQ(n.idOf("nope"), -1);
      CHECK_EQ(n.adj[0].size(), size_t(3)); }  // a-b, a-c, a-d

    // numeric-looking names stay strings and keep insertion order
    { Network n = loadNetwork("test/fixture/numeric.el", false);
      CHECK_EQ(n.names[0], std::string("10"));
      CHECK_EQ(n.names[1], std::string("9"));
      CHECK_EQ(n.idOf("9"), 1); }

    // self-loop is a fatal input error
    { bool threw = false;
      try { loadNetwork("test/fixture/selfloop.el", false); }
      catch (const std::runtime_error&) { threw = true; }
      CHECK(threw); }

    // duplicate edge: last weight wins, edge counted once
    { Network n = loadNetwork("test/fixture/dup.el", true);
      CHECK_EQ(n.m, 1L);
      CHECK_EQ(n.adj[0].size(), size_t(1));
      CHECK_EQ(n.adj[0][0].second, 2.0); }

    CHECK_EQ(choose2(8), 28.0);
    CHECK_EQ(choose2(1), 0.0);

    // n < 2 is a clear error, not a divide by zero
    { bool threw = false;
      try { checkpointNetworkStats(1, 0); }
      catch (const std::runtime_error& e) {
          threw = std::strstr(e.what(), "at least 2 nodes") != nullptr; }
      CHECK(threw); }
}

static void testExternalSort() {
    const size_t kBuf = 256;             // forces many runs
    ExternalSorter s("/tmp/extsort-test", kBuf);
    // 5000 records, shuffled, weight descending should come back ordered.
    for (int i = 0; i < 5000; ++i)
        s.add(std::to_string((i * 7919) % 5000) + " " + std::to_string(i));
    CHECK_EQ(s.size(), 5000L);
    s.finish();
    std::string rec; long count = 0; double prev = 1e300; bool ordered = true;
    while (s.next(rec)) {
        double w = leadingNumber(rec);
        if (w > prev) ordered = false;
        prev = w; ++count;
    }
    CHECK(ordered);
    CHECK_EQ(count, 5000L);
    CHECK_EQ(leadingNumber("0.500000 3 2 1 4"), 0.5);
    CHECK_EQ(leadingNumber("-1 2"), -1.0);
}

int main() {
    testOptions();
    testNetwork();
    testExternalSort();
    return summary();
}
