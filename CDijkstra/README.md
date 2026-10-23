# BLANT-Dijkstra (C++)

BLANT-Dijkstra is a local network alignment algorithm that aligns protein interaction networks from different species and can be used to predict Gene Ontology (GO) annotations. It works in a seed-and-extend manner: it first reads a set of “seeds,” or node pairs, and then grows the alignment outward from them. At each step, a Dijkstra-inspired greedy strategy picks the candidate pair with the best similarity score, while making sure that adding this pair does not violate an edge correctness (EC) threshold, in order to maintain alignment quality. The parameter delta deliberately adds randomness to the choice of candidate, which lets us assess how confident we can be in each aligned pair. The algorithm terminates when there are no more candidate pairs, or when the number of node pairs discarded due to EC violations exceeds a threshold, indicating that no more good pairs exist. Finally, the aligned node pairs (1 to 1) are written to an output file.

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

