#!/usr/bin/env python3
"""
Score community-detection outputs and plot them.

Usage:
  python plot_communities.py --graph edges.txt \
      --run Leiden=leiden.txt --run Infomap=infomap.txt --run SA=sa.txt \
      [--out plot.png] [--include-singletons]

Graph file: one edge per line, "u v" (extra columns such as weights are ignored,
lines starting with # are skipped). Treated as undirected and unweighted;
self-loops and duplicate edges are dropped.

Community file: one community per line, whitespace-separated node IDs.
"""
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_graph(path):
    edges = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in "#%":
                continue
            parts = line.split()
            u, v = int(parts[0]), int(parts[1])
            if u == v:
                continue
            edges.add((min(u, v), max(u, v)))
    return edges


def read_communities(path):
    comms = []
    with open(path) as f:
        for line in f:
            nodes = [int(x) for x in line.split()]
            if nodes:
                comms.append(nodes)
    return comms


def score(edges, comms, include_singletons=False):
    m = len(edges)
    deg = defaultdict(int)
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1

    comm_of = {}
    for i, nodes in enumerate(comms):
        for n in nodes:
            if n in comm_of:
                raise ValueError(f"Node {n} appears in more than one community "
                                 "(overlapping communities are not supported)")
            comm_of[n] = i

    edges_in = defaultdict(int)
    for u, v in edges:
        cu, cv = comm_of.get(u), comm_of.get(v)
        if cu is not None and cu == cv:
            edges_in[cu] += 1

    densities = []
    q = 0.0
    for i, nodes in enumerate(comms):
        n = len(nodes)
        # modularity term: L_c/m - (d_c/2m)^2
        d_c = sum(deg[x] for x in nodes)
        q += edges_in[i] / m - (d_c / (2.0 * m)) ** 2
        # density: internal edges / possible internal edges
        if n < 2:
            if include_singletons:
                densities.append(0.0)
            continue
        densities.append(edges_in[i] / (n * (n - 1) / 2.0))

    avg_density = sum(densities) / len(densities) if densities else 0.0
    return len(comms), avg_density, q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--run", action="append", required=True,
                    help="label=path, repeat for each algorithm")
    ap.add_argument("--out", default="communities_plot.png")
    ap.add_argument("--include-singletons", action="store_true",
                    help="count singleton communities as density 0 (default: skip them)")
    args = ap.parse_args()

    edges = read_graph(args.graph)
    print(f"Graph: {len(edges)} undirected edges")

    results = []
    for spec in args.run:
        label, path = spec.split("=", 1)
        comms = read_communities(path)
        k, dens, q = score(edges, comms, args.include_singletons)
        results.append((label, k, dens, q))
        print(f"{label:>12}: communities={k:5d}  avg edge density={dens:.4f}  Q={q:.4f}")

    fig, ax = plt.subplots(figsize=(7, 5))
    for label, k, dens, _ in results:
        ax.scatter(k, dens, s=80, label=label)
        ax.annotate(label, (k, dens), textcoords="offset points", xytext=(6, 6))
    ax.set_xlabel("Number of communities")
    ax.set_ylabel("Average intra-community edge density")
    ax.set_title("Community detection comparison")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved plot to {args.out}")


if __name__ == "__main__":
    main()
