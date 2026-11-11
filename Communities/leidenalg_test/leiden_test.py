import igraph as ig
import leidenalg as la
from collections import defaultdict
from infomap import Infomap

G = ig.Graph.Read_Edgelist("../../networks/AThaliana.el", directed=False)

# ---------- scoring ----------
def report(G, labels, name):
    deg = G.degree()
    two_m = sum(deg)

    size = defaultdict(int)
    vol = defaultdict(int)
    internal = defaultdict(int)
    cut = defaultdict(int)

    for v, c in enumerate(labels):
        size[c] += 1
        vol[c] += deg[v]
    for u, v in G.get_edgelist():
        cu, cv = labels[u], labels[v]
        if cu == cv:
            internal[cu] += 1
        else:
            cut[cu] += 1
            cut[cv] += 1

    comms = sorted(size)
    cond, dens = [], []
    for c in comms:
        denom = min(vol[c], two_m - vol[c])
        cond.append(cut[c] / denom if denom > 0 else 0.0)
        n = size[c]
        dens.append(internal[c] / (n * (n - 1) / 2) if n > 1 else 0.0)

    Q = G.modularity(labels)  # Newman-Girvan modularity

    '''print(f"=== {name} ({len(comms)} communities) ===")
    print(f"  mean conductance        : {sum(cond)/len(cond):.4f}")
    print(f"  mean intra edge density : {sum(dens)/len(dens):.4f}")
    print(f"  Newman-Girvan Q         : {Q:.4f}\n")'''




if 1:
    leiden = la.find_partition(G, la.ModularityVertexPartition, n_iterations=-1)
    leiden_labels = leiden.membership


    #report(G, leiden_labels, "Leiden (modularity)")

    for line in leiden:
        for node in line:
            print(node, end = " ")
        print()

if 0:

    im = Infomap(seed=123, num_trials=20, silent=True)
    for u, v in G.get_edgelist():
        im.add_link(u, v)
    im.run()
    mods = im.get_modules()
    # vertices with no edges never reach Infomap, so give each its own label
    infomap_labels = [mods.get(n, -1 - n) for n in range(G.vcount())]

    #report(G, infomap_labels, "Infomap")

    communities = defaultdict(list)
    for node_id, module_id in im.get_modules().items():
        communities[module_id].append(node_id)

    for module_id, nodes in sorted(communities.items()):
        print(f"{' '.join(map(str, sorted(nodes)))}")
