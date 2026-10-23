#include "tap.h"
#include "../options.h"
#include "../network.h"
#include "../extsort.h"
#include "../samples.h"
#include "../cluster.h"
#include <cstdio>
#include <sstream>
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
    { Network n = loadNetwork("test/fixture/dup.elw", true);
      CHECK_EQ(n.m, 1L);
      CHECK_EQ(n.adj[0].size(), size_t(1));
      CHECK_EQ(n.adj[0][0].second, 2.0); }

    CHECK_EQ(choose2(8), 28.0);
    CHECK_EQ(choose2(1), 0.0);

    // a one-token line must not silently truncate the rest of the file
    { bool threw = false;
      try { loadNetwork("test/fixture/onetoken.el", false); }
      catch (const std::runtime_error& e) {
          threw = std::strstr(e.what(), "line 2") != nullptr; }
      CHECK(threw); }

    // .el with -w, and .elw without -w, are both refused like the script does
    { bool threw = false;
      try { loadNetwork("test/fixture/fixture.el", true); }
      catch (const std::runtime_error&) { threw = true; }
      CHECK(threw); }
    { bool threw = false;
      try { loadNetwork("test/fixture/dup.elw", false); }
      catch (const std::runtime_error&) { threw = true; }
      CHECK(threw); }

    // the complement round-trips through loadNetwork and is the real complement
    { Network g = loadNetwork("test/fixture/fixture.el", false);
      writeComplement(g, "/tmp/comp-test.el", false);
      Network c = loadNetwork("/tmp/comp-test.el", false);
      CHECK_EQ(c.n, 8L);
      CHECK_EQ(c.m, 15L);
      CHECK_EQ(c.m + g.m, 28L);              // C(8,2)
      std::remove("/tmp/comp-test.el"); }

    // ...and its column count follows the weight flag, so it reloads cleanly
    { Network g = loadNetwork("test/fixture/weighted.elw", true);
      writeComplement(g, "/tmp/comp-test.elw", true);
      Network c = loadNetwork("/tmp/comp-test.elw", true);
      CHECK_EQ(c.n, 4L);
      CHECK_EQ(c.m, 3L);                     // C(4,2) - 3
      std::remove("/tmp/comp-test.elw"); }

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
    // ceil(ED*C(k,2)) - 1, floored at k-1: k=6 has C=15
    CHECK_EQ(minEdgesFor(0.5, 6), 7);   // 7.5 -> ceil 8 -> 7
    CHECK_EQ(minEdgesFor(0.3, 6), 5);   // 4.5 -> ceil 5 -> 4, floored to k-1 = 5
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
    CHECK_EQ(got[0], std::string("100 0 4 1 2 3"));

    // orbit mode keeps every orbit of a passing canon; -h collapses to the best
    std::vector<std::string> orbits, best;
    long no = emitRecords(net, "test/fixture/blant5.orbit.out", m, 'o', 0.5, 5,
                          false, false, [&](const std::string& r){ orbits.push_back(r); });
    long nb = emitRecords(net, "test/fixture/blant5.orbit.out", m, 'o', 0.5, 5,
                          false, true, [&](const std::string& r){ best.push_back(r); });
    CHECK_EQ(no, 16L);
    CHECK_EQ(nb, 8L);

    // BLANT prints %lg counts: fractional under -w, scientific when large
    std::vector<std::string> wrec;
    long nw = emitRecords(net, "test/fixture/blant5.weighted.out", m, 'g', 0.5, 5,
                          false, false, [&](const std::string& r){ wrec.push_back(r); });
    CHECK_EQ(nw, 8L);
    CHECK_EQ(wrec[0], std::string("100.5 0 4 1 2 3"));
    std::vector<std::string> sci;
    long ns = emitRecords(net, "test/fixture/blant5.sci.out", m, 'g', 0.5, 5,
                          false, false, [&](const std::string& r){ sci.push_back(r); });
    CHECK_EQ(ns, 1L);
    CHECK_EQ(leadingNumber(sci[0]), 1000000.0);

    // non-ascending node ids are a corrupt file, not a silent reordering
    bool threw = false;
    try { emitRecords(net, "test/fixture/blant5.bad.out", m, 'g', 0.5, 5, false, false,
                      [](const std::string&){}); }
    catch (const std::runtime_error& e) {
        threw = std::strstr(e.what(), "ascending") != nullptr; }
    CHECK(threw);
}

// Builds a record line "weight node item nbr nbr ..."
static std::string rec(double w, int node, int item, std::vector<int> nbrs) {
    std::ostringstream os; os << w << ' ' << node << ' ' << item;
    for (int v : nbrs) os << ' ' << v;
    return os.str();
}

static void testClusterer() {
    Network net = loadNetwork("test/fixture/fixture.el", false);
    Options opt; opt.overlap = 0.5; opt.minClusArg = 3; opt.maxClus = 1000;
    opt.seed = 1;
    Clusterer c(net, opt, 5, 1.0);
    // Clique {a,b,c,d} = ids 0..3 gets high counts, {e,f,g,h} = 4..7 lower.
    double w[8] = {100, 98, 96, 94, 50, 48, 46, 44};
    std::vector<std::vector<int>> nb = {
        {1,2,3},{0,2,3},{0,1,3},{0,1,2},
        {5,6,7},{4,6,7},{4,5,7},{4,5,6}};
    long idx = 0;
    for (int u = 0; u < 8; ++u)
        c.addRecord(rec(w[u], u, 3, nb[u]), ++idx);
    std::vector<Cluster> got = c.run();
    CHECK(!got.empty());
    const Cluster& best = got[0];
    CHECK(best.nodes.size() >= 3);
    for (size_t i = 1; i < best.nodes.size(); ++i)
        CHECK(best.nodes[i] > best.nodes[i-1]);        // ascending
    CHECK_EQ(best.edges, (long)(best.nodes.size() * (best.nodes.size()-1) / 2));
    CHECK(best.edges > 0);
    CHECK(best.edgeSum > 0);

    // every returned cluster meets the density target
    for (const Cluster& cl : got) {
        CHECK_EQ(cl.edges, (long)(cl.nodes.size()*(cl.nodes.size()-1)/2));
        CHECK_EQ(cl.k, 5);
    }

    // minClus is enforced
    Options small = opt; small.minClusArg = 5;
    Clusterer c2(net, small, 5, 1.0);
    idx = 0;
    for (int u = 0; u < 8; ++u) c2.addRecord(rec(w[u], u, 3, nb[u]), ++idx);
    for (const Cluster& cl : c2.run()) CHECK(cl.nodes.size() >= 5);
}

static Cluster mk(std::vector<int> nodes, int k) {
    Cluster c; c.nodes = nodes; c.k = k;
    size_t s = nodes.size();
    c.edges = (long)(s*(s-1)/2); c.edgeSum = (double)c.edges;
    return c;
}

static void testDedup() {
    // exact duplicates: the second, identical set is dropped
    { std::vector<Cluster> v = { mk({0,1,2},5), mk({0,1,2},5), mk({3,4,5},5) };
      dedupExact(v);
      CHECK_EQ(v.size(), size_t(2)); }

    // a superset arriving later must WIN. Stage 3 sorts by size first, so it
    // removes only exact duplicates and leaves the 3-subset; stage 4 then drops
    // the subset by overlap. The bug this pins: an unsorted stage 3 dropped the
    // 4-clique and kept the 3-clique.
    { std::vector<Cluster> v = { mk({0,1,2},5), mk({0,1,2,3},5) };
      dedupExact(v);
      CHECK_EQ(v.size(), size_t(2));
      CHECK_EQ(v[0].nodes.size(), size_t(4));   // the bigger one leads
      dedupOverlap(v, 0.5);
      CHECK_EQ(v.size(), size_t(1));
      CHECK_EQ(v[0].nodes.size(), size_t(4)); }

    // overlap: 3 of 4 shared is 0.75 > 0.5, dropped
    { std::vector<Cluster> v = { mk({0,1,2,3},5), mk({0,1,2,4},5) };
      dedupOverlap(v, 0.5); CHECK_EQ(v.size(), size_t(1)); }

    // exactly 0.5 overlap is kept: the test is strict >
    { std::vector<Cluster> v = { mk({0,1,2,3},5), mk({0,1,4,5},5) };
      dedupOverlap(v, 0.5); CHECK_EQ(v.size(), size_t(2)); }

    // overlap is measured against the NEW cluster's size only: 3 of 4 shared
    // is 0.75 of the new cluster, so it goes, even though it is only 0.5 of
    // the kept one and a max-size rule would have kept it.
    { std::vector<Cluster> v = { mk({0,1,2,3,4,5},5), mk({0,1,2,6},5) };
      dedupOverlap(v, 0.5); CHECK_EQ(v.size(), size_t(1)); }
}

static void testStartMark() {
    // the script marks i = 1 .. while(i < Slen*(1-overlap)), i.e. ceil(X)-1
    CHECK_EQ(startMarkCount(4, 0.5), 1L);
    CHECK_EQ(startMarkCount(5, 0.5), 2L);
    CHECK_EQ(startMarkCount(3, 0.5), 1L);
    CHECK_EQ(startMarkCount(6, 0.5), 2L);
    CHECK_EQ(startMarkCount(4, 0.0), 3L);
    CHECK_EQ(startMarkCount(4, 1.0), 0L);
}

static void testFormat() {
    Network net = loadNetwork("test/fixture/fixture.el", false);
    Options opt;
    Cluster c = mk({0,1,2,3}, 5);
    CHECK_EQ(formatCluster(c, opt, net),
             std::string("4 nodes, 6/6 edges, 100% density from k5, nodeSet { a b c d }"));
    opt.printMembers = false;
    CHECK_EQ(formatCluster(c, opt, net),
             std::string("4 nodes, 6/6 edges, 100% density from k5"));
    opt.printMembers = true; opt.verbose = true;
    CHECK_EQ(formatCluster(c, opt, net),
             std::string("4 nodes, 6/6 edges, 100% density from k5 from k 5, 1.5 clusterWeight, 6 edgeSum, 1 edgeMean, nodeSet { a b c d }"));
}

int main() {
    testOptions();
    testNetwork();
    testExternalSort();
    testCanonMaps();
    testBlantCommand();
    testEmitRecords();
    testClusterer();
    testDedup();
    testStartMark();
    testFormat();
    return summary();
}
