import clean
import os
import sys
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path + "/snap3.6/")
import snap

def dobfs(graph, x, components, cnum):
    nodes = edges = 0
    queue = [x]
    components[x] = cnum
    while len(queue) > 0:
        src = queue[0]
        nodes += 1
        edges += len(graph[src])
        del queue[0]
        for n in graph[src]:
            if components[n] is None:
                queue.append(n)
                components[n] = cnum
    return nodes, edges//2

def getcomponents(graph, csize):
    n = len(graph)
    components = [None] * n
    cnum = 0
    for i in range(0, n):
        if components[i] is None:
            numnodes, numedges = dobfs(graph, i, components, cnum)
            csize[cnum] = numnodes,numedges
            cnum += 1
    return cnum

def getprops(input_file_name, evalues=None):
    # evalues is the number of eigenvalues to print, when 'None', print all
    output_file_name = "/tmp/" + os.path.basename(input_file_name) + str(os.getpid())

    # node-map
    nodemap = []  # list of lists
    clean.clean(input_file_name, output_file_name, nodemap)  # names inserted into nodemap
    snap_graph = snap.LoadEdgeList(snap.PUNGraph, output_file_name, 0, 1)
    
    edges = 0
    adj_list_graph = [set() for i in range(len(nodemap))]
    for line in open(output_file_name, 'r'):
        line = line.strip()
        a,b = tuple(map(int, line.split()))
        if b not in adj_list_graph[a]: 
            edges += 1
        adj_list_graph[a].add(b)
        adj_list_graph[b].add(a)
    adj_list_graph = [list(x) for x in adj_list_graph]
    
    # delete the temporary file
    os.remove(output_file_name)

    # transitivity (3 * closed_triads/all_triads)
    trans = 0

    nodes = len(nodemap)
    cc_sum = 0  # intermediate sum of node clustering coefficients
    diameter = -1
    ev = []  # list of eigen values
    khop = {}  # k-hops
    degree_dist = {}  # degree hist info

    # components
    component_size = {}  # cnum: (num of nodes, num of edges)
    components = getcomponents(adj_list_graph, component_size)
    csize = list(component_size.values())  # nodes,edges
    csize.sort(reverse=True)

    # for loop for each node
    for nodeid in range(nodes):

        # get node clustering coff
        cc = snap.GetNodeClustCf(snap_graph, nodeid)
        nodemap[nodeid].append(str(cc))
        cc_sum += cc

        # get eccentricity
        ecc = snap.GetNodeEcc(snap_graph, nodeid, False)
        nodemap[nodeid].append(str(ecc))
        diameter = max(diameter, ecc)

        # k-hop (not normalized!)
        nodevec = snap.TIntPr64V()
        snap.GetNodesAtHops(snap_graph, nodeid, nodevec, False)
        for item in nodevec:
            k = item.GetVal1()
            nodes_at_k = item.GetVal2()
            khop[k] = khop.get(k, 0) + nodes_at_k  # increment

    # degree distribution
    degvec = snap.TIntPr64V()
    snap.GetDegCnt(snap_graph, degvec)
    for item in degvec:
        degree_dist[item.GetVal1()] = item.GetVal2()

    # eigen values
    peigv = snap.TFltV()
    if evalues is None:
        evalues = nodes//2
    else:
        evalues = min(evalues//2, nodes//2)  # the SNAP library always computes 2*arg eigenvalues supplied to it. 

    if evalues != 0:
        snap.GetEigVals(snap_graph, evalues, peigv)
        for item in peigv:
            ev.append(item)
        ev.sort(reverse=True)
        for i in range(len(ev)):
            ev[i] = str(format(ev[i], '.4f'))  # 4 digits precision after the decimal point

    # betweenness (node & edge)
    nbw = snap.TIntFlt64H()
    ebw = snap.TIntPrFlt64H()
    snap.GetBetweennessCentr(snap_graph, nbw, ebw, 1.0, False)

    for nodebw in nbw:  # node-betweenness
        nodemap[nodebw].append(str(nbw[nodebw]))

    # OUTPUT!
    print("\n########################################################### Global")
    print("Eigenvalues " + str(len(ev)) + ":", " ".join(ev))
    print("Nodes:", nodes, "Edges:", edges)
    print("Connected components:", components)
    str1 = [str(x[0]) for x in csize]
    str1 = " ".join(str1)
    str2 = [str(x[1]) for x in csize]
    str2 = " ".join(str2)
    print("CC byNode:", str1)
    print("CC byEdge:", str2)
    print("Transitivity (3* triangles/all-triads) [NOT COMPUTED]:", trans)
    print("Avg. clustering coefficient:", cc_sum/nodes)
    print("Diameter:", diameter)

    # K-HOP
    khops = []
    for i in range(1, 1 + max(khop.keys())):
        khops.append(khop.get(i, 0))

    #khops = khop.values()
    khops = list(map(str, khops))
    khops = " ".join(khops)
    print("K-hop distribution:", khops)  # nodes reachable in 1 hop, 2 hops, 3 hops, ....

    # degree distribution
    degree_dist[0] = degree_dist.get(0, 0)
    degrees = []    
    # the 0's
    for i in range(0, 1 + max(degree_dist.keys())):
        degrees.append(degree_dist.get(i, 0))
    
    #degrees = degree_dist.values()
    degrees = list(map(str, degrees))
    degrees = " ".join(degrees)
    print("Degree distribution:", degrees)  # num nodes with degree=0, num nodes with degree=1, num nodes with degree=2, ..... 
    
    #print("NOTE: Betweenness values below (both node & edge) are absolute, and not normailized w.r.t num of nodes in the graph")
    print("nodeName clusCoff eccentricity node_betweenness")
    for nodeid in range(nodes):
        print(" ".join(nodemap[nodeid]))
    print("node1 node2 edge_betweenness")
    for edgebw in ebw:  # edge-betweenness
        print(nodemap[edgebw.GetVal1()][0], nodemap[edgebw.GetVal2()][0], ebw[edgebw])
    print("########################################################### End-of-output\n")


if __name__ == "__main__":    
    error = False
    if len(sys.argv) == 2:
        filename = sys.argv[1].strip()
        getprops(filename)
    elif len(sys.argv) == 4:
        if sys.argv[1] != '-e':
            error = True
        else:
            filename = sys.argv[3].strip()
            evalues = int(sys.argv[2].strip())
            getprops(filename, evalues)
    else:
        error = True

    if error: 
        print("Input error! Usage: python3 main.py [-e eigValsToCompute] inputFile")












