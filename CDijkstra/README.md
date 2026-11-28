# DijSLA (C++)

DijSLA (Dijkstra-inspired Stochastic Local Alignment) is a local network alignment algorithm for protein–protein interaction (PPI) networks of different species. It is designed to use network topology alone, and its alignments can be used to predict Gene Ontology (GO) annotations across species.

DijSLA works in a seed-and-extend manner. It starts from one or more seed node pairs and grows the alignment outward, one node pair at a time. Each candidate pair consists of two unaligned nodes, one from each network, that are neighbors of an already-aligned pair.

Unlike Dijkstra's algorithm, which grows deterministically, DijSLA adds controlled randomness. At each step, let `Smax` be the highest similarity score among the candidates. DijSLA draws a target score `St` uniformly at random from `[Smax - delta, Smax]` and selects the candidate whose score is closest to `St`. The pair is accepted only if the alignment's edge correctness (EC₁ and EC₂) stays at or above the EC threshold; otherwise it is discarded. This lets us sample multiple different high-quality alignments with multiple runs, even from the same seed. Across many runs, node pairs that are aligned more often than expected by chance (high Network Alignment Frequency, NAF) are more likely to be "correct".

The algorithm stops when no candidate pairs remain, or when too many consecutive candidates have been discarded for violating the EC threshold. The resulting one-to-one alignment is written to the output.

---

## How to Run the Code

### 1. Go to the project folder

```bash
cd CDijkstra
```

### 2. Build


Rebuild from scratch:

```bash
make clean && make
```

### 3. Run

The program takes **6 required arguments and 1 optional argument**, in this order:

```bash
./build/cdijkstra <delta> <ECthresh> <simFile> <seedFile> <g1.el> <g2.el> [rngSeed]
```

| # | Argument   | Type   | Description |
|---|------------|--------|-------------|
| 1 | `delta`    | float  | Tolerance parameter that adds randomness during alignment extension (default: `0.05`). |
| 2 | `ECthresh` | double | EC threshold in `[0, 1]` (e.g. `0.1`). |
| 3 | `simFile`  | path   | Node-similarity file between nodes of G1 and G2. |
| 4 | `seedFile` | path   | Seed alignment: node pairs known to match. Accepts any number of seeds (at least one). |
| 5 | `g1.el`    | path   | Edge list for graph 1. |
| 6 | `g2.el`    | path   | Edge list for graph 2. |
| 7 | `rngSeed`  | integer (optional) | Seed for the random number generator. If omitted, a new seed is generated. Every run prints its seed to stderr; passing that value here with the same build and arguments reproduces the run. |


Example:

```bash
./build/cdijkstra 0.05 0.1 data/sim.txt data/seed.txt data/g1.el data/g2.el > alignment.txt 2> log.txt
```

---

## Predefined constants

These are set as macros in `src/main.cpp`. Edit them there, or override them at compile time (for example `make CXXFLAGS+=-DMIN_ITERATIONS=1000`).

| Macro | Default | Description |
|---|---|---|
| `MAX_DISCARD_TIMES` | `10` | Number of consecutive EC-violating candidates after which extension stops. |
| `MIN_ITERATIONS` | `500` | Minimum number of iterations before the `MAX_DISCARD_TIMES` stop rule applies. |

---

## Input file formats

All files are plain text and whitespace-separated (space or tab). Node names are strings.

**Edge list (`g1.el`, `g2.el`)**: one undirected edge per line.
```
A  B
A  C
B  D
```

**Seed file (`seedFile`)**: one pair per line: a node from G1, then its matching node in G2.
```
A  x
B  y
```

**Similarity file (`simFile`)**: one pair per line: a G1 node, a G2 node, and a similarity score.
```
A  x  0.92
A  y  0.10
B  y  0.85
```

Seed and similarity entries whose node names are missing from the matching graph are skipped.

---

## How it works

1. **Parse arguments**: read `delta`, `ECthresh`, and the four file paths.
2. **Map node names**: read each edge list and give every node name an integer index
   (`nameToIndex` / `indexToName`).
3. **Build adjacency matrices**: read each edge list again to fill an
   `N × N` boolean adjacency matrix for each graph.
4. **Read the seed**: convert seed node names into index lists for G1 and G2.
5. **Read similarities**: load a `|V1| × |V2|` similarity matrix (`float`).
6. **Align**: extend the seed alignment using similarity, `delta`, and `ECthresh`.
7. **Write output**: print the aligned pairs to stdout and to the output file.

---

## Output

**stdout / output file**: one aligned pair per line, tab-separated:
```
node1Name	node2Name
```

**stderr**: progress, seed, and timing messages, for example:
```
RNG seed: ...
Seed of graph1: ...
Seed of graph2: ...
[TIME] Reading similarity file: 0.42 seconds
```

