#!/pkg/python/3.7.4/bin/python3
import sys
import math
import heapq
from collections import defaultdict
from graph_helpers import *

def num_graphlets(k):
    if k == 5:
        return 30
    elif k == 4:
        return 9
    else:
        return None

def num_orbits(k):
    if k == 5:
        return 73
    elif k == 4:
        return 15
    else:
        return None

def calc_orbit_counts(k):
    USE_CACHE = True

    if USE_CACHE:
        if k == 5:
            return [1, 2, 2, 2, 3, 4, 3, 3, 4, 3, 4, 4, 4, 4, 3, 4, 6, 5, 4, 5, 6, 6, 4, 4, 5, 5, 8, 4, 6, 6, 7, 5, 6, 6, 6, 5, 6, 7, 7, 5, 7, 7, 7, 6, 5, 5, 6, 8, 8, 6, 6, 8, 6, 9, 6, 6, 4, 6, 6, 8, 9, 6, 6, 8, 8, 6, 7, 7, 8, 5, 6, 6, 4]
        elif k == 4:
            return [1, 2, 2, 2, 3, 4, 3, 3, 4, 3, 4, 4, 4, 4, 3]
        else:
            return None
    else:
        orbit_counts = [None] * num_orbits(k)

        for graphlet_num in range(num_graphlets(k)):
            p = run_orca_raw(k, get_base_graph_path(f'graphlets/graphlet{graphlet_num}'))
            orbit_lines = p.stdout.decode().strip().split('\n')[1:]

            for line in orbit_lines:
                splitted = line.split()
                node_name = splitted[0]
                orbit_num = int(node_name[:-1])
                orbits = splitted[1:]

                # we can't just sum all the orbits, we need to count how many are not zero (because if a node has a degree of 5 we only count that once for "an appearance of orbit 0)
                # we include the orbits effect on itself too
                orbit_count = sum([1 if n != '0' else 0 for n in orbits])

                if orbit_counts[orbit_num] == None:
                    orbit_counts[orbit_num] = orbit_count
                else:
                    assert orbit_counts[orbit_num] == orbit_count

        return orbit_counts

def calc_weights(k):
    orbit_counts = calc_orbit_counts(k)
    weights = [1 - math.log(orbit_count) / math.log(num_orbits(k)) for orbit_count in orbit_counts]
    return weights

class ODVDirectory:
    # file format: every line has node name, followed by orbit counts, separated by spaces
    # NODENAME 23 1 250 37 4 0 ...
    def __init__(self, fname):
        self._directory = dict()

        for line in open(fname, 'r'):
            line_split = line.strip().split()
            node = line_split[0]
            odv_list = [int(s) for s in line_split[1:]]
            odv = ODV(node, odv_list)
            self._directory[node] = odv

    def get_odv(self, node):
        if node in self._directory:
            return self._directory[node]
        else:
            return None

    def get_nodes(self):
        return self._directory.keys()

    def __str__(self):
        return '\n'.join([f'{node}: {odv}' for node, odv in self._directory.items()])


class ODV:
    WEIGHTS = []
    WEIGHT_SUM = 0

    @staticmethod
    def set_weights_vars(k):
        ODV.WEIGHTS = calc_weights(k)
        ODV.WEIGHT_SUM = sum(ODV.WEIGHTS) # 45.08670802954777 <- calculated value from .sim file

    def __init__(self, node, odv_list):
        self._node = node
        self._odv_list = odv_list

    def get_similarity(self, other):
        if len(self._odv_list) == 0 or len(other._odv_list) == 0: # handle the case where the node is not connected to anything in one or both files, causing it to appear with no numbers after it in the .odv file
            return 0

        assert len(self._odv_list) == len(other._odv_list) == len(ODV.WEIGHTS), f'self: {len(self._odv_list)}, other: {len(other._odv_list)}, weights: {len(ODV.WEIGHTS)}, self._node: {self._node}, other._node: {other._node}'
        distance_sum = sum([self._get_single_orbit_similarity(m1, m2, i) for i, (m1, m2) in enumerate(zip(self._odv_list, other._odv_list))])
        weight_sum = ODV.WEIGHT_SUM
        return 1 - distance_sum / weight_sum

    def get_inequal_orbits(self, other):
        assert len(self._odv_list) == len(other._odv_list) == len(ODV.WEIGHTS), f'self: {len(self._odv_list)}, other: {len(other._odv_list)}, weights: {len(ODV.WEIGHTS)}'

        inequal_orbits = []

        for i, (o1, o2) in enumerate(zip(self._odv_list, other._odv_list)):
            if o1 != o2:
                inequal_orbits.append(i)

        return inequal_orbits

    def get_mean_similarity(self, other):
        return mean([self._get_single_orbit_mean_similarity(m1, m2) for m1, m2 in zip(self._odv_list, other._odv_list)])

    def get_odv_val(self, num):
        return self._odv_list[num]

    def __str__(self):
        return ' '.join([str(n) for n in self._odv_list])

    @staticmethod
    def _get_single_orbit_mean_similarity(m1, m2):
        return 1 if m1 == m2 == 0 else min(m1, m2) / max(m1, m2)

    @staticmethod
    def _get_single_orbit_similarity(m1, m2, i):
        # the base of the log doesn't matter
        top_inner = math.log(m1 + 1) - math.log(m2 + 1)
        bot = math.log(max(m1, m2) + 2)
        return ODV.WEIGHTS[i] * abs(top_inner) / bot

def read_in_nodes_wo_deg1(graph_path):
    nodes = read_in_nodes(graph_path)
    adj_set = read_in_adj_set(graph_path)
    new_nodes = [node for node in nodes if len(adj_set[node]) > 1]
    return new_nodes

def get_deg_sim(node1, node2, adj_set1, adj_set2, max_deg1, max_deg2):
    deg1 = len(adj_set1[node1])
    deg2 = len(adj_set2[node2])
    return (deg1 + deg2) / (max_deg1 + max_deg2)

def get_odv_seeds(graph_path1, graph_path2, odv_path1, odv_path2, n, no1=False, alpha=1):
    if no1:
        nodes1 = read_in_nodes_wo_deg1(graph_path1)
        nodes2 = read_in_nodes_wo_deg1(graph_path2)
    else:
        nodes1 = list(read_in_nodes(graph_path1))
        nodes2 = list(read_in_nodes(graph_path2))

    odv_dir1 = ODVDirectory(odv_path1)
    odv_dir2 = ODVDirectory(odv_path2)
    adj_set1 = read_in_adj_set(graph_path1)
    adj_set2 = read_in_adj_set(graph_path2)
    max_deg1 = get_max_deg(adj_set1)
    max_deg2 = get_max_deg(adj_set2)

    top_n = [(-1, '', '')] * n
    heapq.heapify(top_n)
    tot_nodes = len(nodes1) # approximation for less incrementing
    proc_nodes = 0
    percent_printed = 0
    skip = 1

    for node1 in nodes1:
        for i in range(0, len(nodes2), skip):
            node2 = nodes2[i]
            odv1 = odv_dir1.get_odv(node1)
            odv2 = odv_dir2.get_odv(node2)

            if odv1 == None or odv2 == None:
                continue

            odv_sim = odv1.get_similarity(odv2)
            deg_sim = get_deg_sim(node1, node2, adj_set1, adj_set2, max_deg1, max_deg2) # the reason we pass in max is so that we don't have to recalculate it every time we call this
            sim = alpha * odv_sim + (1 - alpha) * deg_sim
            obj = (sim, node1, node2)
            min_top = heapq.heappushpop(top_n, obj)
        
        proc_nodes += 1

        if proc_nodes * 10000 / tot_nodes > percent_printed:
            percent_printed += 1
            # print(f'{proc_nodes} / {tot_nodes}', file=sys.stderr)

    return sorted(top_n, reverse=True)

# function I used to validate the sim function based on Hayes' sim files
def validate_sim_function(gtag1, gtag2):
    FACTOR = 1_000_000

    odv_path1 = get_odv_path(gtag1, 5)
    odv_path2 = get_odv_path(gtag2, 5)
    odv_dir1 = ODVDirectory(odv_path1)
    odv_dir2 = ODVDirectory(odv_path2)
    sim_path = get_fake_seeds_path(f'{gtag1}-{gtag2}', 'sim')

    tot_diff = 0
    tot_pairs = 0
    num_gt10 = 0

    with open(sim_path, 'r') as sim_file:
        for line in sim_file:
            node1, node2, sim = line.strip().split()
            sim = float(sim)
            sim_non_decimal = int(sim * FACTOR) # cuz the sim_path values are rounded to six
            odv1 = odv_dir1.get_odv(node1)
            odv2 = odv_dir2.get_odv(node2)
            my_sim = odv1.get_similarity(odv2)
            my_sim_non_decimal = int(my_sim * FACTOR)
            tot_diff += abs(sim_non_decimal - my_sim_non_decimal)
            tot_pairs += 1

            if tot_pairs % 5000 == 0:
                print(tot_pairs, '/', 1004 ** 2)

            if tot_pairs > 10000:
                break

    avg_diff = (tot_diff / tot_pairs) / FACTOR
    print(f'avg_diff: {avg_diff}')
    print(f'num_gt10: {num_gt10}')

def odv_seeds_to_str(odv_seeds):
    return '\n'.join([f'{node1}\t{node2}\t{score}' for score, node1, node2 in odv_seeds])


if __name__ == '__main__':
    graph_path1 = sys.argv[1]
    graph_path2 = sys.argv[2]
    # generate odvs with orca.g
    odv_path1 = sys.argv[3]
    odv_path2 = sys.argv[4]
    # k is 5 when small enough networks (syeast), 4 for larger networks (iid & temporal)
    k = int(sys.argv[5])
    # n is usually the # of nodes in the smaller of the two networks
    n = int(sys.argv[6])
    ODV.set_weights_vars(k)
    odv_seeds = get_odv_seeds(graph_path1, graph_path2, odv_path1, odv_path2, n)
    print(odv_seeds_to_str(odv_seeds))
