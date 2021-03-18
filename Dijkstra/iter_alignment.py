from collections import defaultdict
import time
from skip_list_iter import *
import random
random.seed(0)
from aligner_helper import *

def fast_align(g1, g2, seed, m, seednum, sims, ec_mode=(0.0, 0.0, 0.0), ed=0.0, sb=0, K=10, delta=0.0,debug=False):
    ## align without skiplist: use array
    alignments = []
    ec1 = ec_mode[0]
    ec2 = ec_mode[1]
    seed1, seed2 = seed

    curralign = Alignment(seed=seed, m=m, seednum=seednum, ec_mode=ec_mode, ed=ed, sb=sb, delta=delta)
    curralign.logfile = g1.name + "_" + g2.name + "_" + str(seednum) + ".log"
    curralign.statsfile = g1.name + "_" + g2.name + "_" + str(seednum) + ".stats"
    for i in range(K):
        start = time.time()
        print("Iteration", i)
        # m is number of edges in seed graphlet
        curralign.g1alignednodes = {seed1}
        curralign.g2alignednodes = {seed2}
        curralign.aligned_pairs = {(seed1, seed2)}
        EA = E1 = E2 = m
        C = get_ec_candidates(g1, g2, seed1, seed2, curralign, ec1, ec2)
        visited = 0
        while len(C) - visited > 0:
            # random pick one candidate from C
            k = random.randint(0, len(C)-visited-1)
            g1node, g2node = C[k]
            # if debug:
            #     print("Pick ", g1node, g2node)

            if g1node in curralign.g1alignednodes or g2node in curralign.g2alignednodes:
                visited += 1
                C[k] = C[-visited]
                continue

            # check sim bound
            if sims[g1node][g2node] < sb:
                visited += 1
                C[k] = C[-visited]
                continue

            # checking ec
            # n1 = num_edges_back_to_subgraph(g1, g1node, curralign.g1alignednodes)
            # n2 = num_edges_back_to_subgraph(g2, g2node, curralign.g2alignednodes)
            # M = num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, curralign.aligned_pairs)
            n1, n2, M = get_n1n2M(g1, g2, g1node, g2node, curralign)
            # print(f"n1={n1}, n2={n2}, M={M}")

            if ((EA + M) / (E1 + n1)) < ec1 or ((EA + M) / (E2 + n2)) < ec2:
                visited += 1
                C[k] = C[-visited]
                continue

            # check edge density
            S_num = len(curralign.aligned_pairs)
            new_ed = (EA + M) / (((S_num + 1) * S_num) / 2)
            if new_ed < ed:
                visited += 1
                C[k] = C[-visited]
                continue

            # adding pair
            # if debug:
            #     print("Adding pair:", (g1node, g2node))
            curralign.g1alignednodes.add(g1node)
            curralign.g2alignednodes.add(g2node)
            curralign.aligned_pairs.add((g1node, g2node))

            print(g1.nodes[g1node], g2.nodes[g2node])

            # update E1, E2, EA
            E1 += n1
            E2 += n2
            EA += M

            # update condidate pairs by exploring with the added pair
            new_C = get_ec_candidates(g1, g2, g1node, g2node, curralign, ec1, ec2)
            if len(new_C) > 0:
                visited += 1
                C[k] = C[-visited]
                C = C[:-visited] + new_C
                visited = 0


        print(time.time() - start)
        print("EA:", EA, "E1:", E1, "E2:", E2)
        alignments.append(curralign.aligned_pairs)
        append_result(g1, g2, curralign)

    return alignments



def fast_align_with_edgefreq(g1, g2, seed, m, seednum, sims, ec_mode=(0.0, 0.0, 0.0), ed=0.0, sb=0, K=10, delta=0.0,debug=False):
    ''' align without skiplist: random pick node to explore; use edge freq to save calculation '''
    alignments = []
    ec1 = ec_mode[0]
    ec2 = ec_mode[1]
    seed1, seed2 = seed

    curralign = Alignment(seed=seed, m=m, seednum=seednum, ec_mode=ec_mode, ed=ed, sb=sb, delta=delta)
    curralign.logfile = g1.name + "_" + g2.name + "_" + str(seednum) + ".log"
    curralign.statsfile = g1.name + "_" + g2.name + "_" + str(seednum) + ".stats"
    for i in range(K):
        start = time.time()
        # m is number of edges in seed graphlet
        curralign.g1alignednodes = {seed1}
        curralign.g2alignednodes = {seed2}
        curralign.aligned_pairs = {(seed1, seed2)}
        EA = E1 = E2 = m
        update_edgefreq(g1, g2, seed1, seed2, curralign, ec1, ec2)

        while len(curralign.edge_freq) > 0:
            # random pick one candidate from C
            g1node, g2node = random.sample(curralign.edge_freq.keys(), 1)[0]

            if g1node in curralign.g1alignednodes or g2node in curralign.g2alignednodes:
                del curralign.edge_freq[(g1node, g2node)]
                print("Removing")
                continue

            # check sim bound
            if sims[g1node][g2node] < sb:
                del curralign.edge_freq[(g1node, g2node)]
                continue

            # n1, n2, M = curralign.edge_freq[(g1node, g2node)]
            n1, n2, M = get_n1n2M(g1, g2, g1node, g2node, curralign)
            if ((EA + M) / (E1 + n1)) < ec1 or ((EA + M) / (E2 + n2)) < ec2:
                del curralign.edge_freq[(g1node, g2node)]
                continue

            # print("Adding pair:", (g1node, g2node))
            # adding pair
            if debug:
                print("Adding pair:", (g1node, g2node))
            curralign.g1alignednodes.add(g1node)
            curralign.g2alignednodes.add(g2node)
            curralign.aligned_pairs.add((g1node, g2node))

            # update E1, E2, EA
            E1 += n1
            E2 += n2
            EA += M

            # update condidate pairs
            update_edgefreq(g1, g2, g1node, g2node, curralign, ec1, ec2)

            del curralign.edge_freq[(g1node, g2node)]


        print("Iteration", i)
        print(time.time() - start)
        alignments.append(curralign.aligned_pairs)
        append_result(g1, g2, curralign)

    return alignments




def iter_align_with_skiplist(g1, g2, seed, m, seednum, sims, ec_mode=(0.0, 0.0, 0.0), ed=0.0, sb=0, K=10, delta=0.0,debug=False):
    # date = str(datetime.datetime.today()[:-10])
    # filename = g1.name+"-"+g2.name+date+".align"
    # to save each alignment
    alignments = []
    ec1 = ec_mode[0]
    ec2 = ec_mode[1]
    seed1, seed2 = seed
    
    C_pairs = get_new_neighbor_pairs(g1, g2, seed1, seed2, {seed1}, {seed2}, sb, sims)
    curralign = Alignment(seed=seed, m=m, seednum=seednum, ec_mode=ec_mode, ed=ed, sb=sb, delta=delta)
    curralign.logfile = g1.name + "_" + g2.name + "_" + str(seednum) + ".log"
    curralign.statsfile = g1.name + "_" + g2.name + "_" + str(seednum) + ".stats"
    for i in range(K):
        start = time.time()
        #m is number of edges in seed graphlet
        curralign.g1alignednodes = {seed1}
        curralign.g2alignednodes = {seed2}
        curralign.aligned_pairs = {(seed1, seed2)}
     
        candidatePairs = set()
        candidatePairs.update(C_pairs)
        curralign.g1candidatenodes = defaultdict(set)
        curralign.g2candidatenodes = defaultdict(set)
        curralign.pq = SkipList()
        curralign.edge_freq = {}
#         update_info(g1, g2, edge_freq, pq, C_pairs, sims, debug)
        E1 = E2 = EA = m
        
        for g1node, g2node in C_pairs:
            n1 = num_edges_back_to_subgraph(g1, g1node, curralign.g1alignednodes)
            n2 = num_edges_back_to_subgraph(g2, g2node, curralign.g2alignednodes)
            M = num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, curralign.aligned_pairs)
            assert(M <= n1 and M <= n2), f"M={M}, n1={n1}, n2={n2}, nodes=({g1node},{g2node})"
            curralign.edge_freq[(g1node, g2node)] = [M, n1, n2]
            curralign.g1candidatenodes[g1node].add(g2node)
            curralign.g2candidatenodes[g2node].add(g1node)
            if debug:
                print((g1node, g2node), " updated in edge_freq ", curralign.edge_freq[(g1node, g2node)])
            val = 2*M/(n1+n2)  # the percentage of the common edges
            curralign.pq.add((val, (g1node, g2node)), debug=debug)
        

        if debug:
            print("aligning inital seeds*************************************************************************")
   
        while len(curralign.pq) > 0:
            # if debug:
            #     curralign.pq.print_list()
            # select the best pair with delta value
            g1node, g2node = best_pair(curralign.pq, delta)
            if debug:
                print("popping out", g1node, g2node)
#                 pq.print_list()
            # look at conditions
            if g1node in curralign.g1alignednodes or g2node in curralign.g2alignednodes:
                continue

            # check sim bound
            if sims[g1node][g2node] < sb:
                continue

            # checking ec
            n1 = num_edges_back_to_subgraph(g1, g1node, curralign.g1alignednodes)
            n2 = num_edges_back_to_subgraph(g2, g2node, curralign.g2alignednodes)
            M = num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, curralign.aligned_pairs)
            
            if ((EA + M)/(E1 + n1)) < ec1 or ((EA + M)/(E2 + n2)) < ec2:
                #print("Fail at ec constraints...continue...")
                continue
            
            # check edge density
            S_num = len(curralign.aligned_pairs)
            new_ed = (EA+M)/(((S_num+1)*S_num)/2)
            if new_ed < ed:
                #print("Fail at ed constraints...continue...")
                continue


            # adding pair
            if debug:            
                print("Adding pair:", (g1node, g2node))
            curralign.g1alignednodes.add(g1node)
            curralign.g2alignednodes.add(g2node)
            curralign.aligned_pairs.add((g1node, g2node))

            # update condidate pairs
            newcandidatePairs = get_new_neighbor_pairs(g1,g2,g1node,g2node, curralign.g1alignednodes, curralign.g2alignednodes, sb, sims)

            # update edge freq
            update_info(g1, g2, curralign, newcandidatePairs, sims, debug)

            # update E1, E2, EA
            E1 += n1
            E2 += n2
            EA += M

        
        print("Iteration",i)
        print(time.time() - start)
        alignments.append(curralign.aligned_pairs)
        append_result(g1, g2, curralign)

    return alignments
