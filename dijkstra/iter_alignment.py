from builder import *
import numpy
from collections import defaultdict
import time, datetime
from skip_list import *
import lzma
import sys
from measures import edgecoverage 
import random


def update_best_pair(pq, yeast_graph, human_graph, yeast_node, human_node, pairs, sims, delta = 0):
##    nonlocal pq
    paired_yeast_nodes = np.fromiter((pair[0] for pair in pairs), dtype=int)
    paired_human_nodes = np.fromiter((pair[1] for pair in pairs), dtype=int)
    yeast_neighbors = np.setdiff1d(
        yeast_graph.get_neighbors(yeast_node), paired_yeast_nodes)
    human_neighbors = np.setdiff1d(
        human_graph.get_neighbors(human_node), paired_human_nodes)
    # do the check here to save on a function call
    # also, saves on unnecessary enqueueing of empty pairs (-1, [])
    # cost: 2 extra comparisons and 1 boolean operation
##    print("yeast", yeast_neighbors)
##    print("human", human_neighbors)
    if yeast_neighbors.size == 0 or human_neighbors.size == 0:
        return
    
    bp_list = sub_best_pair(yeast_neighbors, human_neighbors, sims, delta)
    for (val, new_pairs) in bp_list:
        if val >= 0:
            for pair in new_pairs:
                pq.add((val, pair))
    
                ##pq.insert((val, pair))
    


def update_best_pair(pq, g1, g2, g1node, g2node, pairs, sims, delta = 0):
    paired_g1nodes = np.fromiter((pair[0] for pair in pairs), dtype=int)
    paired_g2nodes = np.fromiter((pair[1] for pair in pairs), dtype=int)
    g1node_neighbors = np.setdiff1d(g1.get_neighbors(g1node), paired_g1nodes)
    g2node_neighbors = np.setdiff1d(g2.get_neighbors(g2node), paired_g2nodes)

    if g1node_neighbors.size == 0 or g2node_neighbors.size == 0:
        return
    
    bp_list = sub_best_pair(g1node_neighbors, g2node_neighbors, sims, delta)
    for (val, new_pairs) in bp_list:
        if val >= 0:
            for pair in new_pairs:
                pq.add((val, pair))
    

def sub_best_pair(yeast_neighbors, human_neighbors, sims, delta):
    """
    sub_best_pair takes takes two sets of yeast and human neighbors
    plus a similarity matrix and returns the highest similarity
    pair of nodes in the two sets of neighbor nodes.
    runtime: O(|Y| * |H|)
    """
    s = sims[np.ix_(yeast_neighbors,human_neighbors)]
    # If the matrix is empty, no pair is found
##    if s.size == 0:
##        return (-1,[])

    #low = 1
    #low = max(s.max() - delta, 0.0)
    low = 0
    answers = []
    while(s.max() >= low):
        y_found, h_found = np.where(s == s.max())
        node_pairs = np.column_stack((yeast_neighbors[y_found], human_neighbors[h_found]))
##    np.random.shuffle(nodes) #do the shuffle when popping instead
        answers.append((s.max(), node_pairs))
        s[y_found, h_found] = -1
    return answers

##    simple dfs is very bad
##    y = random.sample(yeast, 1)[0]
##    h = random.sample(human, 1)[0]
##    return (y,h) + (sims[y][h],)
##
##    (a, b) = random.choice(np.column_stack(found))
##    return (yeast_list[a], human_list[b], s[a][b])
##    list(zip(found[0], found[1])))
##    don't use this for numpy arrays. It runs slower.

def best_pair(pq, delta):
    try:
        pair_list = pq.pop(delta)
    except IndexError:
        raise StopIteration("no more pair values")
    return pair_list[1] 


def fit_in_curr_align(g1, g2, node1, node2, pairs):
    neighbor1 = g1.get_neighbors(node1)
    neighbor2 = g2.get_neighbors(node2)
    flipped = False
    if len(neighbor1) > len(neighbor2):
        g1, g2 = g2, g1
        node1, node2 = node2, node1
        neighbor1, neighbor2 = neighbor2, neighbor1
        flipped = True 
    curr_alignment = {pair[1]:pair[0] for pair in pairs} if flipped else {pair[0]:pair[1] for pair in pairs}

    
    for node in neighbor1:
        if node in curr_alignment and curr_alignment[node] not in neighbor2:
            return False
    return True


def get_new_neighbor_pairs(g1, g2, node1, node2, g1alignednodes, g2alignednodes, sims, sb):
    result = set()
    for i in g1.get_neighbors(node1):
        for j in g2.get_neighbors(node2):
            if i not in g1alignednodes and j not in g2alignednodes:
                if sims[i][j] >= sb:
                    result.add((i,j))
    return result

def get_neighbor_candidate_pairs(g1, g2, node1, node2, g1alignednodes, g2alignednodes, edge_freq, sims):
    new_neighbors = []
    candidate_neighbors = []
    for i in g1.get_neighbors(node1):
        for j in g2.get_neighbors(node2):
            if i not in g1alignednodes and j not in g2alignednodes:
                if (i,j) in edge_freq:
                    candidate_neighbors.append((i,j))            
                else:
                    if sims[i][j] > 0:
                        new_neighbors.append((i,j))

            #if i not in g1alignednodes and j not in g2alignednodes and (i,j) in edge_freq:
            #    candidate_neighbors.append((i,j))            
            #elif i not in g1alignednodes and j not in g2alignednodes:
            #    if sims[i][j] > 0:
            #        neighbors.append((i,j))
    return new_neighbors, candidate_neighbors


def get_neighbor_candidate_pairs2(g1, g2, node1, node2, g1alignednodes, g2alignednodes, g1candidatenodes, g2candidatenodes, edge_freq):
    candidate_neighbors = set()
   
    for g1node in g1.get_neighbors(node1):
        if g1node not in g1alignednodes:
            for g2node in g1candidatenodes[g1node]:
                if (g1node, g2node) in edge_freq:
                    candidate_neighbors.add((g1node,g2node))            

    for g2node in g2.get_neighbors(node2):
        if g2node not in g2alignednodes:
            for g1node in g2candidatenodes[g2node]:
                if (g1node, g2node) in edge_freq:
                    candidate_neighbors.add((g1node,g2node))            

    return candidate_neighbors

def get_aligned_neighbor_pairs(g1,g2, node1,node2, aligned_pairs, trace = False):
    result = []
    
    for i in g1.get_neighbors(node1):
        for j in g2.get_neighbors(node2):
            if (i,j) in aligned_pairs:
                if trace:
                    print(i,j)
                result.append((i,j))
    return result


def strict_align(g1, g2, seed, sims, delta = 0):
    used_node1 = set()
    used_node2 = set()
    stack = []
    pairs = set()
    # initialize
    #print(g1.nodes)
    for seed1, seed2 in seed:
        #pairs.add((seed1, seed2, sims[seed1][seed2]))
        #print(g1.nodes[seed1])
        #print(seed1)
        pairs.add((seed1, seed2, sims[seed1][seed2]))
        used_node1.add(seed1)
        used_node2.add(seed2)
        stack += get_neighbor_pairs(g1,g2,seed1,seed2,sims) 

    # while we still have still nodes to expand
    while stack:
        node1, node2 = stack.pop()
        if node1 not in used_node1 and node2 not in used_node2:
            used_node1.add(node1)
            used_node2.add(node2)
            if fit_in_curr_align(g1,g2,node1,node2,pairs):
                pairs.add((node1,node2,sims[node1][node2]))
                #pairs.add((seed1, seed2, sims[g1.indexes[seed1]][g2.indexes[seed2]]))
                stack += get_neighbor_pairs(g1,g2,node1,node2,sims)
    return (used_node1, used_node2, pairs)


def stop_align2(g1, g2, seed, sims, ec_mode, delta = 0):
    g1alignednodes = set()
    g2alignednodes = set()
    aligned_pairs = set()
    pq = SkipList()
    stack = []
    
    for seed1, seed2 in seed:
        aligned_pairs.add((seed1, seed2, sims[seed1][seed2]))
        g1alignednodes.add(seed1)
        g2alignednodes.add(seed2)
        stack += get_neighbor_pairs(g1,g2,seed1,seed2,sims) 
        update_best_pair(pq, g1, g2, seed1, seed2, aligned_pairs, sims, delta)


    while len(g1alignednodes) < len(g1):
        try:
            curr_pair = best_pair(pq, delta)
            g1node = curr_pair[0]
            g2node = curr_pair[1]
            if g1node in g1alignednodes or g2node in g2alignednodes:
                continue
            if fit_in_curr_align(g1, g2, g1node, g2node, aligned_pairs):
                update_best_pair(pq, g1, g2, g1node, g2node, aligned_pairs, sims, delta)
                aligned_pairs.add((g1node, g2node, sims[g1node][g2node]))
                g1alignednodes.add(g1node)
                g2alignednodes.add(g2node)
        except(StopIteration): 
            break

    return (g1alignednodes, g2alignednodes, aligned_pairs)

def num_edges_back_to_subgraph(graph, node, aligned_nodes, trace=False):
    edges = 0
   
    for neighbor_node in aligned_nodes:
        if graph.has_edge(node, neighbor_node):
            if trace:
                print("Nnode in graph : " + str(neighbor_node)) 
            edges += 1

    '''
    for neighbor_node in graph.get_neighbors(node):
        if neighbor_node in aligned_nodes:
            edges += 1
    '''
    return edges


def num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, aligned_pairs):
    edgepairs = 0
    #print("printing edges for M")
    for n1, n2 in aligned_pairs:
        if g1.has_edge(g1node, n1) and g2.has_edge(g2node, n2):
            #print(n1, n2)
            edgepairs += 1

    """
    for n1 in g1.get_neighbors(g1node):
        for n2 in g2.get_neighbors(g2node):
            if (n1, n2) in aligned_pairs:
                edgepairs += 1
    """
    return edgepairs


def update_edge_freq(g1, g2, candidatePairs, g1alignednodes, g2alignednodes, edge_freq):
    for g1node, g2node in candidatePairs:
        n1 = num_edges_back_to_subgraph(g1, g1node, g1alignednodes)   
        n2 = num_edges_back_to_subgraph(g2, g2node, g2alignednodes)   
        M = num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, aligned_pairs) 
        edge_freq[(g1node, g2node)] = [M, n1, n2]

    return edge_freq

def iter_align(g1, g2, seed, sims, ec_mode, ed, m, sb=0, K=10, debug=False):
    # Need to change the output filename
    filename = g1.name+"-"+g2.name+str(random.randint(0,10000))+".align"
    # to save each alignment
    alignments = []
    ec1 = ec_mode[0]
    ec2 = ec_mode[1]
    seed1, seed2 = seed

    C_pairs = get_new_neighbor_pairs(g1,g2,seed1,seed2,{seed1}, {seed2},sims, sb)
    for i in range(K):
        #m is number of edges in seed graphlet
        g1alignednodes = {seed1}
        g2alignednodes = {seed2}
        aligned_pairs = {(seed1, seed2)}
    
        candidatePairs = set()
        candidatePairs.update(C_pairs) 
        E1 = E2 = EA = m
        #edge_freq = {}

        if debug:
            print("aligning inital seeds*************************************************************************")
   
        while len(candidatePairs) > 0:
            print(len(candidatePairs))
            # not doing combination checking now
            # use default python set.pop(); may need to use a more random pop() later
            cand_list = list(candidatePairs)
            g1node, g2node = cand_list.pop(random.randrange(len(candidatePairs)))
            candidatePairs.remove((g1node, g2node))
            
            # look at conditions
            if g1node in g1alignednodes or g2node in g2alignednodes:
                continue

            # check sim bound
            if sims[g1node][g2node] < sb:
                continue

            # checking ec
            n1 = num_edges_back_to_subgraph(g1, g1node, g1alignednodes)   
            n2 = num_edges_back_to_subgraph(g2, g2node, g2alignednodes)   
            M = num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, aligned_pairs) 
            
            if ((EA + M)/(E1 + n1)) < ec1 or ((EA + M)/(E2 + n2)) < ec2:
                continue
            
            # check edge density
            S_num = len(aligned_pairs)
            new_ed = (EA+M)/(((S_num+1)*S_num)/2)
            if new_ed < ed:
                continue
            
            #print("Adding pair:", (g1node, g2node))

            # adding pair
            g1alignednodes.add(g1node)
            g2alignednodes.add(g2node)
            aligned_pairs.add((g1node, g2node))

            # update condidate pairs
            candidatePairs.update(get_new_neighbor_pairs(g1,g2,g1node,g2node, g1alignednodes, g2alignednodes,sims, sb)) 

 
            # update edge freq
            #edge_freq[(g1node, g2node)] = [M, n1, n2]
            #g1candidatenodes[g1node].add(g2node)
            #g2candidatenodes[g2node].add(g1node)

            # update E1, E2, EA
            E1 += n1
            E2 += n2
            EA += M
        
        print("Iteration",i)
        alignments.append(aligned_pairs)
        write_result(filename, aligned_pairs, g1, g2)
            
    return alignments
            



def output(k, E1, E2, EA, seed, runtime, seednum, size):
    print("seednum: " + str(seednum) + " k:" + str(k) +  " size:" + str(size) + " E1:" + str(E1) + " E2:" + str(E2) + " EA:" + str(EA) + " time:" + str(runtime) + " seed: " + str(seed))



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

#def s3score(g1, g2, pairs, subgraph):
#    aligned_edges = len(subgraph)
#    u2in = unaligned_edges_g2_in(g1, g2, pairs, subgraph)
#    denom = aligned_edges + (g1.num_edges() / 2) + len(u2in) 
#    return aligned_edges/denom


def ec1score(E1, EA):
    return EA/E1

def ec2score(E2, EA):
    return EA/E2

def s3score(E1,E2, EA):
    return EA/(E1+E1-EA)    

def write_result(filename, pairs, graph1, graph2):
    with open(filename, 'a+')as f:
        #f.write(str(len(d)) + ' ' + str(coverage(yeast_graph, human_graph,d)) + '\n')
        for x in pairs:
            print(str(graph1.nodes[x[0]]) + ' ' + str(graph2.nodes[x[1]]))
            f.write(str(graph1.nodes[x[0]]) + ' ' + str(graph2.nodes[x[1]]) + '\n')

        f.write('\n')


def to_name(pairs, yd, hd):
    return [(yd[yeast_graph], hd[human_graph]) for (yeast_graph, human_graph, sims) in pairs]


#import datetime
def log_file_name(start = 'bionet_yeast_human', ext = '.txt'):
    dtime = datetime.datetime.now()
    return start + '_' +dtime.strftime("%Y-%m-%d_%H-%M-%S") + ext
