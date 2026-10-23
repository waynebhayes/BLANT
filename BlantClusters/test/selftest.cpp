#include "tap.h"
#include "../options.h"
#include "../network.h"
#include "../extsort.h"
#include "../samples.h"
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

static void testCanonMaps() {
    CanonMaps m;
    CHECK(m.load(5, "test/fixture"));
    CHECK_EQ(m.canonEdges.size(), size_t(34));
    CHECK_EQ(m.canonEdges[0], 0);      // empty graphlet
    CHECK_EQ(m.canonEdges[1], 1);      // one edge
    CHECK_EQ(m.canonEdges[3], 3);
    CHECK_EQ(m.orbit2canon[0], 0);     // orbit_map line 2 -> canon 0
    CHECK_EQ(m.orbit2canon[4], 2);     // orbit 4 belongs to canon 2
    CHECK_EQ(m.orbit2canon.size(), size_t(90));
    CHECK_EQ(m.orbit2canon[9], 4);
    CHECK(!m.load(7, "test/fixture")); // missing files

    // minimum edge count: int(ED * C(k,2)) - 1, floored at k-1
    CHECK_EQ(minEdgesFor(0.5, 5), 4);   // int(5.0)-1 = 4
    CHECK_EQ(minEdgesFor(1.0, 5), 9);   // int(10.0)-1 = 9
    CHECK_EQ(minEdgesFor(0.1, 5), 4);   // int(1.0)-1 = 0, floored to k-1 = 4
    CHECK_EQ(minEdgesFor(0.6, 5), 5);   // int(6.0)-1 = 5
}

static void testBlantCommand() {
    Options o;
    o.blantCmd = "./blant"; o.blantExe = "./blant";
    CHECK_EQ(blantCommand(o, 5, "net.el"),
             std::string("./blant -p0.9L -q -k5 -sEBE -mcg net.el"));
    Options w; w.blantCmd = "./blant"; w.blantExe = "./blant";
    w.weighted = true; w.seedArg = "-r 42";
    CHECK_EQ(blantCommand(w, 4, "net.elw"),
             std::string("./blant -p0.9L -q -w -r 42 -k4 -sEBE -mcg net.elw"));
    Options v; v.blantCmd = "./blant"; v.blantExe = "./blant";
    v.verbose = true; v.quiet = false; v.complement = true; v.turan = true;
    v.seedArg = "-r 7";
    CHECK_EQ(blantCommand(v, 6, "c.el", "12"),
             std::string("./blant -p0.9L -A 12 -C -r 7 -k6 -sEBE -mcg c.el"));
    CHECK_EQ(blantCommand(v, 6, "c.el", ""),
             std::string("./blant -p0.9L -C -r 7 -k6 -sEBE -mcg c.el"));
}

static void testEmitRecords() {
    Network net = loadNetwork("test/fixture/fixture.el", false);
    CanonMaps m; CHECK(m.load(5, "test/fixture"));
    std::vector<std::string> got;
    long n = emitRecords(net, "test/fixture/blant5.out", m, 'g', 0.5, 5,
                         false, false, [&](const std::string& r){ got.push_back(r); });
    CHECK_EQ(n, 8L);
    CHECK_EQ(got.size(), size_t(8));
    CHECK_EQ(got[0], std::string("40 0 4 1 2 3"));

    // orbit mode keeps every orbit of a passing canon; -h collapses to the best
    std::vector<std::string> orbits, best;
    long no = emitRecords(net, "test/fixture/blant5.orbit.out", m, 'o', 0.5, 5,
                          false, false, [&](const std::string& r){ orbits.push_back(r); });
    long nb = emitRecords(net, "test/fixture/blant5.orbit.out", m, 'o', 0.5, 5,
                          false, true, [&](const std::string& r){ best.push_back(r); });
    CHECK_EQ(no, 16L);
    CHECK_EQ(nb, 8L);

    // non-ascending node ids are a corrupt file, not a silent reordering
    bool threw = false;
    try { emitRecords(net, "test/fixture/blant5.bad.out", m, 'g', 0.5, 5, false, false,
                      [](const std::string&){}); }
    catch (const std::runtime_error& e) {
        threw = std::strstr(e.what(), "ascending") != nullptr; }
    CHECK(threw);
}

int main() {
    testOptions();
    testNetwork();
    testExternalSort();
    testCanonMaps();
    testBlantCommand();
    testEmitRecords();
    return summary();
}
