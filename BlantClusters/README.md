# blant-clusters (C++)

A C++17 port of `scripts/blant-clusters.sh`. It is an interface-compatible
replacement: same command line and output format. Results differ from the
script because of three intentional bug fixes and a different RNG (see below).
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

# Overview

 ```
   BlantClusters/
     blant-clusters.cpp   main: tmp dir + signal cleanup, stage sequencing
     options.{h,cpp}      the script's option switch, verbatim rules
     network.{h,cpp}      .el/.elw loading, complement writer
     samples.{h,cpp}      canon maps, `-mc` parsing, BLANT invocation
     extsort.{h,cpp}      bounded-memory external merge sort
     cluster.{h,cpp}      BFS community growth, pruning, dedup, output
     test/                selftest (unit), check (invariant), e2e.sh (golden)
     makefile
 ```

## Dependencies
It invokes blant itself; BLANT is not reimplemented. 
No libwayne dependency

# Pipeline

1. Network. Read the edge list, map names to dense int ids, keep adjacency as (neighbour, weight).
2. BLANT. Runs once per k, stdout to the temp dir. Canon maps come from canon_maps/canon_list$k.txt and orbit_map$k.txt (canon is the line number minus 2, edge count is field 3, orbits map to canons through orbit_map)
3. Stage 1. Stream the BLANT output, keep canons with `minEdges = max(ceil(ED·C(k,2))−1, k−1)`, accumulate counts and neighbour lists one node at a time, and emit records of the form weight node item nbr….
4. Sort. Variable-length records, byte budget `-B` (512 MB default), sorted runs spilled to temp files, then a k-way heap merge.
5. Stage 2. Greedy BFS over the graphlet-neighbour graph with the script's `count[v]/count[u] ≥ 0.9√ED` admission rule, its `max(|S|, n/100)` miss budget, the 3σ and mode-peel prune, then a peel until the density target is reached.
6. Stages 3 and 4. Exact duplicates are removed with bigger sets first (dropping superset copy). Then overlap dedup.

# Differences from the script

The port fixes three script bugs on purpose, and each fix changes which clusters are found:
1. `AppendNeighbors(u, orign)` (line 348): the typo turned off the `line[v] > line[origin]` restriction for every node after the first.
2. `-h` (lines 213, 220): it set every item to 0 and then skipped item 0, so the script printed nothing.
3. `-M12` (line 110): the `sed` strips `-m`, not `-M`.
4. stage-1 sort ties break on record text
5. the final peel recomputes the weakest node after each removal (the script uses one in-degree snapshot)
6. the BFS shuffle uses `std::mt19937` instead of gawk's `randsort`.

The comparison is on CElegans (`networks/CElegans.el`) with one shared BLANT run (`./blant -n 200000 -r 42 -k5 -sEBE -mcg`), given to both programs with `-D`, seeds 1–3. "fixed script" is the script with only the `orign` typo corrected. Cells show clusters / distinct nodes covered:

| ED  | script            | fixed script      | C++               |
|-----|-------------------|-------------------|-------------------|
| 1   | 4 / 16            | 1 / 4             | 1 / 4 (identical bytes) |
| 0.9 | 6 / 28–29         | 4 / 19            | 4 / 18            |
| 0.7 | 30–35 / 123–140   | 22–26 / 92–103    | 23–25 / 92–97     |
| 0.5 | 242–266 / 556–580 | 224–235 / 530–558 | 225–247 / 522–542 |

Most of the gap from the original script is the `orign` fix. Every cluster from both programs meets its density target, but the stricter rules find fewer of them. 

# Performance

BLANT runs for all k in parallel, like the script (two k on CElegans: 159% CPU). The clustering stages run one (ED, k) at a time and take 0.03 s on CElegans against the script's 0.6 s. The weighted adjacency index is built once per run.

# Tests

1. `test/selftest` covers option parsing, network edge cases, the external sort with input larger than its buffer, canon maps, BLANT command assembly, record emission, cluster growth, dedup, start-mark arithmetic and output formatting.
2. `test/check` re-derives density, edge mean and overlap from the network alone and knows nothing about how the clusters were found.
3. `test/e2e.sh` checks the golden file, same-seed byte determinism.
4. e2e runs a fake BLANT without `-D` (command line, quoting with a space in the path, redirect, non-zero exit), checks a weighted run, and feeds the checker a cluster line over 64 KB.
7. `regression-tests/blantClusters/test.sh` builds and tests the port under `regression-test-all.sh`.

# Known limit

The default `-p0.9L` is rejected by the current `./blant` (`unknown _stopMode 1`) for the script and the port alike; use `-n N`.
