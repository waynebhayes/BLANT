#!/bin/python3
from index_helpers import *
from graph_helpers import *
from general_helpers import *

class PatchedIndexEntry:
    def __init__(self, entry1, entry2, adj_set):
        # left is min right is max
        if entry1.get_key() <= entry2.get_key():
            left_entry = entry1
            right_entry = entry2
        else:
            left_entry = entry2
            right_entry = entry1

        self._left_connectors = left_entry.get_matching_list(right_entry) # returns matching nodes in left's array poses
        self._right_connectors = right_entry.get_matching_list(left_entry) # returns matching nodes in right's array poses
        self._num_matching = len(self._left_connectors)

        # node order in _patched_node_arr [nonmatching nodes of left_entry in order] + [matching nodes in left_entry order] + [nonmatching nodes of right_entry in order]
        # this is the canonical ordering, assuming both graphlets are multiplicity 1. that means, when checking if the seed is a perfect ortholog match, we will make sure that all nodes in one patched entry's _patched_node_arr are orthologs with the node in the same position in the other patched entry's _patched_node_arr
        self._patched_node_arr = []
        self._matching_poses = [] # stores (g1_pos, g2_pos) of matching nodes
        matching_left_poses = set()
        matching_right_poses = set()

        for i, node in enumerate(left_entry.get_node_arr()):
            if not right_entry.node_in(node):
                self._patched_node_arr.append(node)

        for i, node in enumerate(left_entry.get_node_arr()):
            if right_entry.node_in(node):
                self._patched_node_arr.append(node)

                assert right_entry.index_of(node) != -1

                self._matching_poses.append((i, right_entry.index_of(node)))
                matching_left_poses.add(i)
                matching_right_poses.add(right_entry.index_of(node))

        for i, node in enumerate(right_entry.get_node_arr()):
            if not left_entry.node_in(node):
                self._patched_node_arr.append(node)

        # create extra_edges
        self._extra_edges = dict()

        non_matching_left_poses = [i for i in range(len(left_entry)) if i not in matching_left_poses]
        non_matching_right_poses = [i for i in range(len(right_entry)) if i not in matching_right_poses]
        
        assert len(right_entry) == len(left_entry)
        assert len(non_matching_left_poses) == len(non_matching_right_poses) == len(right_entry) - self._num_matching

        for left_pos in non_matching_left_poses:
            for right_pos in non_matching_right_poses:
                if left_entry.get_node_arr()[left_pos] in adj_set[right_entry.get_node_arr()[right_pos]]:
                    if left_pos not in self._extra_edges:
                        self._extra_edges[left_pos] = []

                    self._extra_edges[left_pos].append(right_pos)

        extra_edges_str = '-'.join([f'{key}:' + f'{",".join(str(pos) for pos in value)}' for key, value in sorted(list(self._extra_edges.items()))])
        matching_poses_str = ','.join([f'{g1_pos}:{g2_pos}' for g1_pos, g2_pos in self._matching_poses])
        self._patch_id = f'{left_entry.get_key()}-{right_entry.get_key()};{matching_poses_str};{extra_edges_str}'

    def get_key(self):
        return self._patch_id

    def get_node_arr(self):
        return self._patched_node_arr
                    
    def __str__(self):
        return self._patch_id


def get_matching_poses_list(entry_list, prox, target_num_matching):
    matching_poses_list = []
    total_to_check = len(entry_list)
    num_checked = 0

    for i in range(0, len(entry_list)):
        matching_poses_list.append([])
        outer_entry = entry_list[i]

        for j in range(i + prox, i + prox + 1):
            if j >= len(entry_list):
                continue

            inner_entry = entry_list[j]
            curr_num_matching = outer_entry.get_num_matching(inner_entry)

            if target_num_matching < 0:
                do_append = curr_num_matching >= -1 * target_num_matching
            else:
                do_append = curr_num_matching == target_num_matching

            if do_append:
                matching_poses_list[i].append(j)

            num_checked += 1

    return matching_poses_list

def get_patched_index(k, index_path, graph_path, prox=6, target_num_matching=6):
    entry_list = read_in_entry_list(index_path, k)
    matching_poses_list = get_matching_poses_list(entry_list, prox, target_num_matching)
    adj_set = read_in_adj_set(graph_path)
    patched_index = Index()

    for i, matching_list in enumerate(matching_poses_list):
        for j in matching_list:
            patched_entry = PatchedIndexEntry(entry_list[i], entry_list[j], adj_set)
            patched_index.add_entry(patched_entry)

    return patched_index

if __name__ == '__main__':
    patched_index = get_patched_index(8, get_index_path('mouse'), get_graph_path('mouse'), 6, 6)
    print(len(patched_index))
    assert_with_prints(len(patched_index), 6833, 'len(patched_index)')
