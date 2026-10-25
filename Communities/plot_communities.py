#!/usr/bin/env python3
"""
Plot every community from a single run: x = community size, y = edge density.

Usage:
  python plot_communities.py --graph edges.txt --run Leiden=leiden.txt \
      [--out plot.png] [--logx] [--include-singletons]

--run accepts "label=path" or just "path".

Graph file: one edge per line, "u v" (extra columns ignored, lines starting
with # or % skipped). Treated as undirected and unweighted; self-loops and
duplicate edges are dropped.

Community file: one community per line, whitespace-separated node IDs.
"""
import argparse
import os
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


def per_community(edges, comms, include_singletons=False):
    """Return (sizes, densities, Q). Singletons are skipped unless requested."""
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

    sizes, densities = [], []
    q = 0.0
    for i, nodes in enumerate(comms):
        n = len(nodes)
        d_c = sum(deg[x] for x in nodes)
        q += edges_in[i] / m - (d_c / (2.0 * m)) ** 2
        if n < 2:
            if include_singletons:
                sizes.append(n)
                densities.append(0.0)
            continue
        sizes.append(n)
        densities.append(edges_in[i] / (n * (n - 1) / 2.0))
    return sizes, densities, q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph", required=True)
    ap.add_argument("--run", required=True, help="label=path or path (single run)")
    ap.add_argument("--out", default="communities_plot.png")
    ap.add_argument("--logx", action="store_true", help="log-scale x axis")
    ap.add_argument("--include-singletons", action="store_true",
                    help="plot singleton communities with density 0 (default: skip)")
    args = ap.parse_args()

    if os.path.exists(args.run) or "=" not in args.run:
        label, path = os.path.splitext(os.path.basename(args.run))[0], args.run
    else:
        label, path = args.run.split("=", 1)

    edges = read_graph(args.graph)
    comms = read_communities(path)
    sizes, densities, q = per_community(edges, comms, args.include_singletons)

    mean_d = sum(densities) / len(densities) if densities else 0.0
    print(f"Graph: {len(edges)} undirected edges")
    print(f"{label}: communities={len(comms)}  plotted={len(sizes)}  "
          f"mean density={mean_d:.4f}  Q={q:.4f}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(sizes, densities, s=18, alpha=0.6, edgecolors="none")
    ax.axhline(mean_d, color="gray", linestyle="--", linewidth=1,
               label=f"mean density = {mean_d:.3f}")
    if args.logx:
        ax.set_xscale("log")
    ax.set_xlabel("Community size (nodes)")
    ax.set_ylabel("Edge density")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"{label}: {len(sizes)} communities")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved plot to {args.out}")


if __name__ == "__main__":
    main()
