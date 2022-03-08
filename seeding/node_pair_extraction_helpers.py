#!/bin/python3
from collections import defaultdict

def extract_node_pairs(all_seeds_list):
    node_pair_voting = create_node_pair_voting(all_seeds_list)
    node_favorite_pairs = create_node_favorite_pairs(node_pair_voting)
    output_pairs = create_output_pairs(node_favorite_pairs)
    return output_pairs

def aug(node, n):
    return f'{n}_{node}'

def deaug(auged_node):
    return '_'.join(auged_node.split('_')[1:])

def print_output_pairs(output_pairs):
    print('\n'.join([f'{deaug(node1)} {deaug(node2)}' for node1, node2 in output_pairs]))

def create_node_pair_voting(all_seeds_list):
    def add_to_voting(node1, node2):
        if node1 not in node_pair_voting:
            node_pair_voting[node1] = defaultdict(int)

        if node2 not in node_pair_voting:
            node_pair_voting[node2] = defaultdict(int)

        node_pair_voting[node1][node2] += 1
        node_pair_voting[node2][node1] += 1

    node_pair_voting = dict()

    for graphlet_id, s1_index, s2_index in all_seeds_list:
        for s1_node, s2_node in zip(s1_index, s2_index):
            add_to_voting(aug(s1_node, 1), aug(s2_node, 2))

    return node_pair_voting

def create_node_favorite_pairs(node_pair_voting):
    node_favorite_pairs = defaultdict(set)

    for base, votes in node_pair_voting.items():
        max_count = max([count for count in votes.values()])

        for node, count in votes.items():
            if count == max_count:
                node_favorite_pairs[base].add(node)

    return node_favorite_pairs

def create_output_pairs(node_favorite_pairs):
    output_pairs = set()

    for node, favorites in node_favorite_pairs.items():
        for fav in favorites:
            if node == fav:
                exit('node equals fav')

            if node < fav: # only process in one direction to avoid duplicates
                if node in node_favorite_pairs[fav]:
                    output_pairs.add((deaug(node), deaug(fav)))

    return output_pairs
