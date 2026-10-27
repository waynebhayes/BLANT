# blant-clusters (C++)

A C++17 port of `scripts/blant-clusters.sh`: same command line, same stdout
format. It still runs `./blant` for the graphlet sampling; only the awk/sort
pipeline around it is reimplemented.

**Status:** `scripts/blant-clusters.sh` stays the authoritative version, and
`scripts/do-blant-runs.sh`, `make-clusters.sh` and `blant-clusters-predict.sh`
keep calling it until this port has been compared on more networks. Make an
algorithm change in both places until then.

## Build and test

    make            # builds ./blant-clusters (needs only a C++17 compiler)
    make test       # unit tests, invariant checker, end-to-end tests

`regression-tests/blantClusters/test.sh` runs `make all test` as part of
`regression-test-all.sh`. Run the binary from the BLANT root, like the script:
it reads `canon_maps/` and runs BLANT relative to the current directory.

## Intentional differences from the script

These are script bugs, fixed on purpose. Each one changes which clusters are
found, so the output does not match the script's exactly:

1. `AppendNeighbors(u, orign)` (script line 348): the typo made `line[origin]`
   empty for every node after the first, which silently turned off the
   `line[v] > line[origin]` restriction. The port applies it as written.
2. `-h` (script lines 213, 220): it set every item to 0 and then skipped item 0,
   so the script printed nothing. The port keeps the best item per node.
3. `-M12` (script line 110): the `sed` strips `-m`, not `-M`, so the value was
   never parsed. The port reads it.

Other differences: ties in the stage-1 sort are broken by record text, and the
BFS neighbour order uses `std::mt19937`, not gawk's `randsort`. So for the
same `-r` seed the two programs explore in a different order. Cluster counts
and sizes come out in the same range (see below), but individual clusters can
differ.

## Comparison with the script

CElegans (`networks/CElegans.el`), one shared BLANT run
(`./blant -n 200000 -r 42 -k5 -sEBE -mcg`), both programs given it with `-D`,
seeds `-r 1..3`. "fixed script" is the script with only the `orign` typo
corrected. Cells show clusters / distinct nodes covered:

| ED  | script             | fixed script    | C++             |
|-----|--------------------|-----------------|-----------------|
| 1   | 4 / 16             | 1 / 4           | 1 / 4 (same bytes as fixed) |
| 0.9 | 6 / 28–29          | 4 / 19          | 4 / 18          |
| 0.7 | 30–35 / 123–140    | 22–26 / 92–103  | 23–25 / 92–97   |
| 0.5 | 242–266 / 556–580  | 224–235 / 530–558 | 225–247 / 522–542 |

The ED 0.9 gap is one cluster, `{2 21 22 25 26}` (9/10 edges). The C++ finds
it with 3 of 10 seeds and the script with all 10. Whether node 26 is admitted
depends only on whether the BFS dequeues it before node 22.

## Performance

- BLANT runs for all k start at once and are waited for together, like the
  script, so wall-clock time is the slowest k, not the sum. Peak memory is the
  sum of the concurrent BLANT runs, as in the script.
- The clustering stages run one (ED, k) pair at a time. On the CElegans run
  above they take 0.03 s against the script's 0.6 s.
- `-B` (default 512 MB) bounds the in-memory part of the stage-1 sort. All
  spill files are open at once during the merge, so keep the run count
  (intermediate bytes / `-B`) under the open-file limit.

## Known limits

- The default `-p0.9L` precision is passed to BLANT unchanged. The current
  `./blant` rejects it (`RunBlantFromGraph: unknown _stopMode 1`), for the
  script as well as this port; pass `-n N` instead.
