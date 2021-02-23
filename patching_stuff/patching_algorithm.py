import sys
import json
import re

# functional constants
MISSING_ALLOWED = 0

# command line input
def get_graph_path_for_species(species):
    return f'/home/sana/Jurisica/IID/networks/IID{species}.el'

species1 = sys.argv[1]
species2 = sys.argv[2]
s1_index_file = open(sys.argv[3], 'r')
s2_index_file = open(sys.argv[4], 'r')
NUM_MATCHING_NODES = int(sys.argv[5]) if len(sys.argv) >= 6 else 6 # negative for <=, positive for ==
PATCH_PROX_INC = int(sys.argv[6]) if len(sys.argv) >= 7 else 1
SEED_PROX_INC = int(sys.argv[7]) if len(sys.argv) >= 8 else 1
DEBUG = bool(eval(sys.argv[8])) if len(sys.argv) >= 9 else True

# setup
def debug_print(*args, **kwargs):
    if DEBUG:
        print(*args, file = sys.stderr, **kwargs)

s1_graph_file = open(get_graph_path_for_species(species1), 'r')
s2_graph_file = open(get_graph_path_for_species(species2), 'r')

# main code
class Index:
    def __init__(self, index_str):
        splitted_str = index_str.split(' ')
        self._graphlet_id = int(splitted_str[0])
        self._node_arr = splitted_str[1:]
        self._node_set = set(splitted_str[1:])

    def get_graphlet_id(self):
        return self._graphlet_id

    def get_node_arr(self):
        return self._node_arr

    def node_in(self, node):
        return node in self._node_set

    def get_num_matching(self, other_index):
        num_matching = 0

        for node in self._node_arr:
            if node in other_index._node_set:
                num_matching += 1

        return num_matching

    def get_matching_list(self, other_index):
        matching_list = []

        for i, node in enumerate(self._node_arr):
            if node in other_index._node_set:
                matching_list.append(i)

        return matching_list

    def __str__(self):
        return ' '.join(self._node_arr)

    def __len__(self):
        return len(self._node_arr)


class PatchedIndex:
    def __init__(self, index1, index2, adj_set):
        # left is min right is max
        if index1.get_graphlet_id() <= index2.get_graphlet_id():
            left_index = index1
            right_index = index2
        else:
            left_index = index2
            right_index = index1

        self._left_connectors = left_index.get_matching_list(right_index) # returns matching nodes in left's array poses
        self._right_connectors = right_index.get_matching_list(left_index) # returns matching nodes in right's array poses
        self._num_matching = len(self._left_connectors)

        # node order is [eight nodes of left_index in order] + [nodes of right_index not in left_index, in order]
        self._patched_node_arr = []
        self._matching_poses = [] # stores g1 poses of matching nodes

        for i, node in enumerate(left_index.get_node_arr()):
            if not right_index.node_in(node):
                self._patched_node_arr.append(node)

        for i, node in enumerate(left_index.get_node_arr()):
            if right_index.node_in(node):
                self._patched_node_arr.append(node)
                self._matching_poses.append(i)

        for node in right_index.get_node_arr():
            if not left_index.node_in(node):
                self._patched_node_arr.append(node)

        # create extra_edges
        self._extra_edges = dict()
        
        for i in range(len(left_index) - self._num_matching):
            for j in range(len(left_index), len(self._patched_node_arr)):
                if self._patched_node_arr[j] in adj_set[self._patched_node_arr[i]]:
                    if i not in self._extra_edges:
                        self._extra_edges[i] = []

                    self._extra_edges[i].append(j)

        extra_edges_str = ';'.join([f'{key}:' + f'{",".join(str(pos) for pos in value)}' for key, value in sorted(list(self._extra_edges.items()))])
        matching_poses_str = ','.join([str(pos) for pos in self._matching_poses])
        self._key = f'{left_index.get_graphlet_id()};{right_index.get_graphlet_id()};{matching_poses_str};{extra_edges_str}'

    def get_node_arr(self):
        return self._patched_node_arr
                    
    def __str__(self):
        return self._key


def read_graph_file(graph_file):
    adj_set = dict()

    for edge_str in graph_file:
        node1, node2 = re.split('[\s\t]', edge_str.strip())

        if node1 not in adj_set:
            adj_set[node1] = set()

        if node2 not in adj_set:
            adj_set[node2] = set()

        adj_set[node1].add(node2)
        adj_set[node2].add(node1)

    return adj_set

def read_index_file(index_file):
    index_list = []

    for curr_index_str in index_file:
        curr_index = Index(curr_index_str.strip())
        index_list.append(curr_index)

    return index_list

def get_matching_poses_adj_list(index_list):
    PRINT = False # this is different than DEBUG because it controls whether or not we should print this section of the code, given DEBUG == True
    matching_poses_adj_list = []
    total_to_check = len(index_list) * PATCH_PROX_INC
    num_checked = 0
    percent_printed = 0

    for i in range(0, len(index_list)):
        matching_poses_adj_list.append([])
        outer_index = index_list[i]

        for j in range(i + 1, i + 1 + PATCH_PROX_INC):
            if j >= len(index_list):
                continue

            inner_index = index_list[j]

            if NUM_MATCHING_NODES < 0:
                do_append = outer_index.get_num_matching(inner_index) >= -1 * NUM_MATCHING_NODES
            else:
                do_append = outer_index.get_num_matching(inner_index) == NUM_MATCHING_NODES

            if do_append:
                matching_poses_adj_list[i].append(j)

            num_checked += 1

        if PRINT:
            if num_checked * 100 / total_to_check > percent_printed:
                debug_print(f'{percent_printed}% done')
                percent_printed += 1

    num_total_matches = sum(len(matching_list) for matching_list in matching_poses_adj_list)
    debug_print(f'for {NUM_MATCHING_NODES} matching nodes, there are {num_total_matches} out of {num_checked} matches, or {num_total_matches * 100 / num_checked}%')
    return matching_poses_adj_list

def get_patched_indexes(matching_poses_adj_list, index_list, adj_set):
    patched_indexes = dict()

    for i, matching_list in enumerate(matching_poses_adj_list):
        for j in matching_list:
            patched_index = PatchedIndex(index_list[i], index_list[j], adj_set)
            patched_index_key = str(patched_index)
            
            if patched_index_key not in patched_indexes:
                patched_indexes[patched_index_key] = []

            patched_indexes[patched_index_key].append(patched_index)

    return patched_indexes

def get_s1_to_s2():
    ortho_file = open('/home/wayne/src/bionets/SANA/Jurisica/IID/Orthologs.Uniprot.tsv', 'r')
    SPECIES_TO_INDEX = dict()
    species_line = ortho_file.readline().strip()
    species_order = re.split('[\s\t]+', species_line)

    for i, species in enumerate(species_order):
        SPECIES_TO_INDEX[species] = i

    s1_to_s2 = dict()
    s1_pos = SPECIES_TO_INDEX[species1]
    s2_pos = SPECIES_TO_INDEX[species2]

    for line in ortho_file:
        line_split = line.strip().split()

        if line_split[s1_pos] == species1: # first line
            assert line_split[s2_pos] == species2
        else: # other lines
            s1_node = line_split[s1_pos]
            s2_node = line_split[s2_pos]

            if s1_node != '0' and s2_node != '0':
                s1_to_s2[s1_node] = s2_node

    debug_print(f'there are {len(s1_to_s2.values())} orthologs and {len(set(s1_to_s2.values()))} unique orthologs')
    return s1_to_s2

def get_orthologs_list(s1_indexes, s2_indexes, s1_to_s2, num_seeds):
    PRINT = False
    orthologs_list = []
    pairs_processed = 0
    total_pairs_to_process = len(s1_indexes) * SEED_PROX_INC ** 2
    percent_printed = 0

    for patched_id, s1_patched_indexes in s1_indexes.items():
        if patched_id in s2_indexes:
            for i in range(0, min(len(s1_patched_indexes), SEED_PROX_INC)):
                s1_index = s1_patched_indexes[i].get_node_arr()

                for j in range(0, min(len(s2_indexes[patched_id]), SEED_PROX_INC)):
                    s2_index = s2_indexes[patched_id][j].get_node_arr()
                    assert(len(s1_index) == len(s2_index) and len(s1_index) > 8), f's1_index:{s1_index}, s2_index:{s2_index}'
                    missing_nodes = 0

                    for m in range(len(s2_index)):
                        s1_node = s1_index[m]
                        s2_node = s2_index[m]

                        if s1_node not in s1_to_s2 or s1_to_s2[s1_node] != s2_node:
                            missing_nodes += 1

                            if missing_nodes > MISSING_ALLOWED:
                                break

                    if missing_nodes <= MISSING_ALLOWED:
                        orthologs_list.append((patched_id, s1_index, s2_index))

                    pairs_processed += 1

                # print
                if PRINT:
                    if pairs_processed / total_pairs_to_process * 100 > percent_printed:
                        debug_print(f'{percent_printed}% done', file=sys.stderr)
                        percent_printed += 1

    print(f'on settings NUM_MATCHING_NODES={NUM_MATCHING_NODES} PATCH_PROX_INC={PATCH_PROX_INC} SEED_PROX_INC={SEED_PROX_INC}, there are {len(orthologs_list)} {MISSING_ALLOWED}|miss orthologs out of {pairs_processed} processed patched pairs, representing {len(orthologs_list) * 100 / pairs_processed}%')
    return orthologs_list

def main():
    # read in index files
    s1_adj_set = read_graph_file(s1_graph_file)
    s1_index_list = read_index_file(s1_index_file)
    s2_adj_set = read_graph_file(s2_graph_file)
    s2_index_list = read_index_file(s2_index_file)

    # find pairs with enough nodes matching
    s1_matching_poses_adj_list = get_matching_poses_adj_list(s1_index_list)
    s2_matching_poses_adj_list = get_matching_poses_adj_list(s2_index_list)

    # patch indexes
    s1_patched_indexes = get_patched_indexes(s1_matching_poses_adj_list, s1_index_list, s1_adj_set)
    s2_patched_indexes = get_patched_indexes(s2_matching_poses_adj_list, s2_index_list, s2_adj_set)

    # find num seeds
    num_seeds = 0

    for key in s1_patched_indexes:
        if key in s2_patched_indexes:
            num_seeds += len(s1_patched_indexes[key]) * len(s2_patched_indexes[key])

    debug_print(f'there are {num_seeds} total possible seeds (patched index pairs)')

    # calculate ortholog percentage
    s1_to_s2 = get_s1_to_s2()
    orthologs_list = get_orthologs_list(s1_patched_indexes, s2_patched_indexes, s1_to_s2, num_seeds)


if __name__ == '__main__':
    main()
