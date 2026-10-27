# blant-clusters (C++)

A C++17 port of `scripts/blant-clusters.sh`. It is an interface-compatible
replacement: same command line and output format. Results differ from the
script in two option fixes and in the BFS's random order (see below).
It still runs `./blant` for the graphlet sampling; only the awk/sort pipeline
around it is reimplemented.

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

## Differences from the script

Two script bugs are fixed on purpose:

1. `-h` (script lines 213, 220): it set every item to 0 and then skipped item 0,
   so the script printed nothing. The port keeps the best item per node.
2. `-M12` (script line 110): the `sed` strips `-m`, not `-M`, so the value was
   never parsed. The port reads it.

The script's actual algorithm is kept, including these behaviours a
reader might take for bugs:

- `AppendNeighbors(u, orign)` (script line 348) is a typo, but it means the
  `line[v] > line[origin]` restriction binds only the origin's own neighbours.
  Applying it to every node finds fewer clusters at the same density (at ED 1
  on CElegans: 1 cluster of 4 nodes instead of 4 clusters covering 16), so the
  port keeps the script's behaviour. Whether to tighten it is an algorithm
  question for the maintainer.
- The final peel removes nodes in the order of one in-degree snapshot
  (smallest first, ties by node name), not by recomputing the weakest node.
- Stage-1 records carry node names, and ties in the sort are broken by the
  whole record, descending, like `sort -gr` under `LC_ALL=C`.

What still differs: the BFS neighbour order uses `std::mt19937`, not gawk's
`randsort`. The script's final dedup orders same-size clusters by whichever
node gawk's hash order lists first, which cannot be reproduced. Ties in the
script's sort also follow your locale unless `LC_ALL=C` is set. Its output
labels every cluster with the k of the last cluster it read; the port prints each
cluster's own k.

## Comparison with the script

CElegans (`networks/CElegans.el`), one shared BLANT run
(`./blant -n 200000 -r 42 -k5 -sEBE -mcg`) given to both with `-D`, seeds
`-r 1..6`, script under `LC_ALL=C`. Mean clusters / distinct nodes covered:

| ED  | script        | C++           |
|-----|---------------|---------------|
| 1   | 4.0 / 16.0    | 4.0 / 16.0    |
| 0.9 | 6.0 / 29.0    | 6.3 / 29.8    |
| 0.7 | 32.8 / 130.8  | 32.8 / 127.2  |
| 0.5 | 247.7 / 564.8 | 241.3 / 557.8 |
| 0.3 | 107.3 / 444.0 | 107.3 / 432.5 |

Up to ED 0.5 the gaps are within one standard error over the six seeds (node
coverage varies by 7–14 between seeds). At ED 0.3 the C++ covers 2.6% fewer
nodes, about two standard errors. The likely cause is the dedup order above,
but this has not been confirmed.

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
