from pathlib import Path
from collections import defaultdict
from skip_list import *
import uuid

class Alignment:
    def __init__(self, seed, m, ec_mode=(0.0, 0.0, 0.0), ed=0.0, sb=0.0, alpha=0.0, delta=0.0, seednum=0, outputdir="",
                 timestop=-1, alignstop=1000):
        self.g1alignednodes = set()
        self.g2alignednodes = set()
        self.aligned_pairs = set()
        self.edge_freq = {}
        self.g1candidatenodes = defaultdict(set)
        self.g2candidatenodes = defaultdict(set)
        self.pq = SkipList()

        self.numaligns = 0
        self.maxalignsize = 0
        self.alignstop = alignstop
        self.ec1 = ec_mode[0]
        self.ec2 = ec_mode[1]
        self.ed = ed
        self.sb = sb
        self.E1 = m
        self.E2 = m
        self.EA = m
        self.g1seedstr = ""
        self.g2seedstr = ""
        self.k = 0
        self.seed = seed
        self.alpha = alpha
        self.delta = delta
        self.seednum = seednum
        self.currtime = 0
        self.recdepth = 0
        self.timestop = timestop
        self.g1name = ""
        self.g1name = ""
        self.outputdir = outputdir
        self.logfile = ""
        self.statsfile = ""


def best_pair(pq, delta, fast=True):
    '''
    Pop out a pair in the skiplist within range [max-delta, max];
    fast=True means pop with random value but the first node
    '''
    try:
        pair_list = pq.pop(delta, fast)
        # print(pair_list)
    except IndexError:
        raise StopIteration("no more pair values")
    return pair_list[1]

############ funtions to get the neighbors ##############

def get_new_neighbor_pairs(g1, g2, node1, node2, g1alignednodes, g2alignednodes, simbound, sims):
    '''
    Return possible neighbor pairs from given (node1, node2) pair.
    '''
    result = set()
    for i in g1.get_neighbors(node1):
        for j in g2.get_neighbors(node2):
            if i not in g1alignednodes and j not in g2alignednodes:
                if sims[i][j] >= simbound:
                    result.add((i, j))
    return result


def get_candidate_pairs(g1, g2, node1, node2, g1alignednodes, g2alignednodes, g1candidatenodes, g2candidatenodes,
                        edge_freq, simbound, sims):
    '''
    return the new candidate pairs and the existed candidate pairs
    '''
    newcp = set()
    existingcp = set()

    for i in g1.get_neighbors(node1):
        for j in g2.get_neighbors(node2):
            g1nodealigned = i in g1alignednodes
            g2nodealigned = j in g2alignednodes
            if not g1nodealigned and not g2nodealigned:
                if sims[i][j] >= simbound:
                    newcp.add((i, j))
            elif not g1nodealigned:
                for g2node in g1candidatenodes[i]:
                    if (i, g2node) in edge_freq:
                        existingcp.add((i, g2node))
            elif not g2nodealigned:
                for g1node in g2candidatenodes[j]:
                    if (g1node, j) in edge_freq:
                        existingcp.add((g1node, j))

    return newcp, existingcp


def get_neighbor_candidate_pairs(g1, g2, node1, node2, g1alignednodes, g2alignednodes, g1candidatenodes,
                                 g2candidatenodes, edge_freq, simbound, sims):
    '''
    return the exisiting neighbor candidatePairs
    '''
    candidate_neighbors = set()

    for g1node in g1.get_neighbors(node1):
        if g1node not in g1alignednodes:
            for g2node in g1candidatenodes[g1node]:
                if (g1node, g2node) in edge_freq:
                    if sims[g1node][g2node] >= simbound:
                        candidate_neighbors.add((g1node, g2node))

    for g2node in g2.get_neighbors(node2):
        if g2node not in g2alignednodes:
            for g1node in g2candidatenodes[g2node]:
                if (g1node, g2node) in edge_freq:
                    if sims[g1node][g2node] >= simbound:
                        candidate_neighbors.add((g1node, g2node))

    return candidate_neighbors


def get_aligned_neighbor_pairs(g1, g2, node1, node2, aligned_pairs, trace=False):
    result = []
    for i in g1.get_neighbors(node1):
        for j in g2.get_neighbors(node2):
            if (i, j) in aligned_pairs:
                if trace:
                    print(i, j)
                result.append((i, j))
    return result

def get_ec_candidates(g1, g2, g1node, g2node, curralign, ec1, ec2):
    candidate_neighbors = []
    EA, E1, E2 = curralign.EA, curralign.E1, curralign.E2
    for node1 in g1.get_neighbors(g1node):
        for node2 in g2.get_neighbors(g2node):
            if node1 in curralign.g1alignednodes or node2 in curralign.g2alignednodes:
                continue
            n1, n2, M = get_n1n2M(g1, g2, node1, node2, curralign)
            if ((EA + M) / (E1 + n1)) >= ec1 and ((EA + M) / (E2 + n2)) >= ec2:
                candidate_neighbors.append((node1, node2))
    return candidate_neighbors


############ Functions for calculate n1, n2, M (to satisfy ec constraints)

def num_edges_back_to_subgraph(graph, node, aligned_nodes, trace=False):
    edges = 0
    for neighbor_node in aligned_nodes:
        if graph.has_edge(node, neighbor_node):
            if trace:
                print(str(node), " in graph : " + str(neighbor_node))
            edges += 1
    return edges


def num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, aligned_pairs):
    edgepairs = 0
    for n1, n2 in aligned_pairs:
        if g1.has_edge(g1node, n1) and g2.has_edge(g2node, n2):
            edgepairs += 1
    return edgepairs

def get_n1n2M(g1, g2, g1node, g2node, curralign):
    n1 = num_edges_back_to_subgraph(g1, g1node, curralign.g1alignednodes)
    n2 = num_edges_back_to_subgraph(g2, g2node, curralign.g2alignednodes)
    M = num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, curralign.aligned_pairs)
    assert (
                M <= n1 and M <= n2), f"M={M}, n1={n1}, n2={n2}, nodes=({g1node},{g2node}), curralign=({curralign.aligned_pairs})"
    return n1, n2, M


##### helper functions for updating infomation related to curralign after cp insertion

def update_info(g1, g2, curralign, candidatePairs, sims, debug):
    for g1node, g2node in candidatePairs:
        if g1node in curralign.g1alignednodes or g2node in curralign.g2alignednodes:
            if debug:
                print("updating edge_freq: ", (g1node, g2node), " already aligned")
            # delete this pair from edge_freq
            # TODO make the val calculation in a seperate function
            val = 2 * curralign.edge_freq[(g1node, g2node)][0] / (
                        curralign.edge_freq[(g1node, g2node)][1] + curralign.edge_freq[(g1node, g2node)][2])
            curralign.pq.remove_by_name((val, (g1node, g2node)))
            del curralign.edge_freq[(g1node, g2node)]
            continue
        if sims[g1node][g2node] < curralign.sb:
            continue
        if (g1node, g2node) in curralign.edge_freq:
            val = 2 * curralign.edge_freq[(g1node, g2node)][0] / (
                    curralign.edge_freq[(g1node, g2node)][1] + curralign.edge_freq[(g1node, g2node)][2])
            curralign.pq.remove_by_name(val, (g1node, g2node))
            del curralign.edge_freq[(g1node, g2node)]
        n1 = num_edges_back_to_subgraph(g1, g1node, curralign.g1alignednodes)
        n2 = num_edges_back_to_subgraph(g2, g2node, curralign.g2alignednodes)
        M = num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, curralign.aligned_pairs)
        assert (M <= n1 and M <= n2 and (n1 > 0 or n2 > 0)), f"M={M}, n1={n1}, n2={n2}, nodes=({g1node},{g2node})"
        curralign.edge_freq[(g1node, g2node)] = [M, n1, n2]
        curralign.g1candidatenodes[g1node].add(g2node)
        curralign.g2candidatenodes[g2node].add(g1node)
        pair = (g1node, g2node)
        if debug:
            print(pair, " updated in edge_freq ", curralign.edge_freq[pair])
        # val = curralign.edge_freq[pair][0]
        val = 2 * M / (n1 + n2)  # the percentage of the common edges
        curralign.pq.add((val, pair), debug=debug)


def update_edgefreq(g1, g2, g1node, g2node, curralign, ec1, ec2):
    EA, E1, E2 = curralign.EA, curralign.E1, curralign.E2
    for node1 in g1.get_neighbors(g1node):
        for node2 in g2.get_neighbors(g2node):
            if node1 in curralign.g1alignednodes or node2 in curralign.g2alignednodes:
                continue
            n1, n2, M = get_n1n2M(g1, g2, node1, node2, curralign)
            if ((EA + M) / (E1 + n1)) < ec1 or ((EA + M) / (E2 + n2)) < ec2:
                continue
            curralign.edge_freq[(node1, node2)] = (n1, n2, M)


def update_skip_list(g1, g2, curralign, candidatePairs, sims, debug):
    for g1node, g2node in candidatePairs:
        if g1node in curralign.g1alignednodes or g2node in curralign.g2alignednodes:
            if debug:
                print("updating edge_freq: ", (g1node, g2node), " already aligned")
            #delete this pair from edge_freq?
            #del curralign.edge_freq[(g1node, g2node)]
            continue
        n1 = num_edges_back_to_subgraph(g1, g1node, curralign.g1alignednodes)   
        n2 = num_edges_back_to_subgraph(g2, g2node, curralign.g2alignednodes)   
        M = num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, curralign.aligned_pairs)            
        assert(M <= n1 and M <= n2), f"M={M}, n1={n1}, n2={n2}, nodes=({g1node},{g2node})"
        curralign.edge_freq[(g1node, g2node)] = [M, n1, n2]
        curralign.g1candidatenodes[g1node].add(g2node)
        curralign.g2candidatenodes[g2node].add(g1node)
        pair = (g1node, g2node)
        if debug:
            print(pair, " updated in edge_freq ", curralign.edge_freq[pair])
        # val = curralign.edge_freq[pair][0]
        val = 2 * M / (n1 + n2)  # the percentage of the common edges
        curralign.pq.add((val,pair), debug=debug)


####### functions to present the results

def induced_subgraph(graph1, graph2, aligned_pairs):
    result = []
    for p in aligned_pairs:
        result.extend(
            [((p[0], q[0]), (p[1], q[1]))
                 for q in aligned_pairs
                    if graph1.has_edge(p[0],q[0])
                        and graph2.has_edge(p[1],q[1])])
    return result

def induced_graph1(graph1, aligned_pairs):
    result = []
    for p in aligned_pairs:
        result.extend(
            [((p[0], q[0]), (p[1], q[1]))
                 for q in aligned_pairs
                    if graph1.has_edge(p[0],q[0])])
    return result

def induced_graph2(graph2, aligned_pairs):
    result = []
    for p in aligned_pairs:
        result.extend(
            [((p[0], q[0]), (p[1], q[1]))
                 for q in aligned_pairs
                    if graph2.has_edge(p[1],q[1])])
    return result

def coverage(yeast, human, subgraph):
    y = len(subgraph) / (yeast.num_edges() / 2) * 100.0
    h = len(subgraph) / (human.num_edges() / 2) * 100.0
    return (y,h)

def unaligned_edges_g2_in(graph1, graph2, aligned_pairs, subgraph):
    uedges = []
    while aligned_pairs:
        p = aligned_pairs.pop()
        uedges.extend( [ (p[1], q[1]) for q in aligned_pairs if graph2.has_edge(p[1], q[1]) and ((p[0], q[0]), (p[1], q[1])) not in subgraph ]   )
    return uedges

def ec1score(E1, EA):
    return EA/E1

def ec2score(E2, EA):
    return EA/E2

def s3score(E1,E2, EA):
    return EA/(E1+E1-EA)


####### Functions to write/print results/log

def writelog(g1,g2, curralign):
    g1edges = induced_graph1(g1, curralign.aligned_pairs)
    g2edges = induced_graph2(g2, curralign.aligned_pairs)
    eaedges = induced_subgraph(g1, g2, curralign.aligned_pairs)
    print("ec1score: ", ec1score(len(g1edges) / 2, len(eaedges) / 2))
    print("ec2score: ", ec2score(len(g2edges) / 2, len(eaedges) / 2))
    print("s3score: ", s3score(len(g1edges) / 2, len(g2edges) / 2, len(eaedges) / 2))


def writestats(curralign):
    f = open(curralign.statsfile, "w")
    f.write("seed: " + curralign.g1seedstr + "\n")
    f.write("k: " + curralign.k + "\n")
    f.write("EC1 bound: " + curralign.ec1 + "\n")
    f.write("EC2 bound: " + curralign.ec2 + "\n")
    f.write("ED bound: " + curralign.ed + "\n")
    f.write("timestop: " + curralign.timestop + "\n")


def printoutput(g1, g2, curralign):
    size = len(curralign.aligned_pairs)
    hours, rem = divmod(curralign.currtime, 3600)
    minutes, seconds = divmod(rem, 60)
    runtime = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)

    result = ("seednum: " + str(curralign.seednum) + " k:" + str(curralign.k) + " size:" + str(size) + " E1:" + str(
        curralign.E1) + " E2:" + str(curralign.E2) + " EA:" + str(curralign.EA) + " time:" + str(
        runtime) + " seed: " + curralign.g1seedstr)

    if curralign.outputdir == "":
        output_dir = "seed" + str(curralign.seednum)
    else:
        output_dir = curralign.outputdir + "/seed" + str(curralign.seednum)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_file = output_dir + "/" + curralign.logfile

    with open(output_file, "a") as f:
        f.write(result + "\n")
        for n1, n2 in curralign.aligned_pairs:
            f.write(str(g1.nodes[n1]) + " " + str(g2.nodes[n2]) + "\n")


def write_result(g1, g2, curralign):
    uuidstr = str(uuid.uuid4())
    uid = uuidstr[:13]
    fname = g1.name + "--" + g2.name + "--" + str(curralign.delta) + "--" + str(curralign.k) + "--" + uid + ".dijkstra"
    if curralign.outputdir == "":
        output_dir = "seed" + str(curralign.seednum)
    else:
        output_dir = curralign.outputdir + "/seed" + str(curralign.seednum)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_file = output_dir + "/" + fname
    with open(output_file, 'w+')as f:
        for x in curralign.aligned_pairs:
            print(str(g1.nodes[x[0]]) + ' ' + str(g2.nodes[x[1]]))
            f.write(str(g1.nodes[x[0]]) + ' ' + str(g2.nodes[x[1]]) + '\n')


def append_result(g1, g2, curralign):
    # uuidstr = str(uuid.uuid4())
    # uid = uuidstr[:13]
    fname = g1.name + "--" + g2.name + "--" + str(curralign.k) + "--seed" + str(curralign.seednum) + ".dijkstra"

    if curralign.outputdir == "":
        output_dir = "seed" + str(curralign.seednum)
    else:
        output_dir = curralign.outputdir + "/seed" + str(curralign.seednum)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    output_file = output_dir + "/" + fname
    with open(output_file, 'a')as f:
        for x in curralign.aligned_pairs:
            # print(str(g1.nodes[x[0]]) + ' ' + str(g2.nodes[x[1]]))
            f.write(str(g1.nodes[x[0]]) + ' ' + str(g2.nodes[x[1]]) + '\n')
        f.write('\n')






