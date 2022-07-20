import re
from general_helpers import *
from species_helpers import *

# does not include syeast
def all_species():
    return ['cat', 'chicken', 'cow', 'dog', 'duck', 'guinea_pig', 'horse', 'human', 'mouse', 'pig', 'rabbit', 'rat', 'sheep', 'turkey']

def get_graph_path(species):
    return f'{get_seeding_dir()}/networks/{species_to_full_name(species)}.el'

def read_in_adj_set(graph_path):
    with open(graph_path, 'r') as graph_file:
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

def read_in_el(graph_fname):
    el = []
    graph_file = open(graph_fname, 'r')

    for line in graph_file:
        node1, node2 = re.split('[\s\t]', line.strip())
        el.append((node1, node2))

    graph_file.close()
    return el

def read_in_nodes(graph_path):
    nodes = set()
    graph_file = open(graph_path, 'r')

    for edge_str in graph_file:
        node1, node2 = re.split('[\s\t]', edge_str.strip())
        nodes.add(node1)
        nodes.add(node2)

    graph_file.close()
    return nodes
