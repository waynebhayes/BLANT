#!/pkg/python/3.7.4/bin/python3
from collections import defaultdict

# extracts node pairs from many2many alignments (.aln files) 
def extract_node_pairs_from_m2m(m2m_pairs):
    node_pair_voting = create_node_pair_voting(m2m_pairs)
    node_favorite_pairs = create_node_favorite_pairs(node_pair_voting)
    output_pairs = create_output_pairs(node_favorite_pairs)
    return output_pairs

def aug(node, n):
    return f'{n}_{node}'

def deaug(auged_node):
    return '_'.join(auged_node.split('_')[1:])

def seeds_to_m2m(seeds):
    # has to be list, not set, because we want duplicates (they count towards the vote)
    m2m_pairs = list()

    for graphlet_id, s1_index, s2_index in seeds:
        for s1_node, s2_node in zip(s1_index, s2_index):
            m2m_pairs.append((s1_node, s2_node))

    return m2m_pairs

def extract_node_pairs(all_seeds_list):
    m2m_pairs = seeds_to_m2m(all_seeds_list)
    pairs = extract_node_pairs_from_m2m(m2m_pairs)
    return pairs

'''def get_rid_of_deg1_pairs(pairs, gtag1, gtag2):
    graph_path1 = get_graph_path(gtag1)
    nodes1 = read_in_nodes(graph_path1)
    adj_set1 = read_in_adj_set(graph_path1)
    graph_path2 = get_graph_path(gtag2)
    nodes2 = read_in_nodes(graph_path2)
    adj_set2 = read_in_adj_set(graph_path2)
    deg1_nodes1 = [node for node in nodes1 if len(adj_set1[node]) == 1]
    deg1_nodes2 = [node for node in nodes2 if len(adj_set2[node]) == 1]
    new_pairs = []

    for node1, node2 in pairs:
        if node1 not in deg1_nodes1 and node2 not in deg1_nodes2:
            new_pairs.append((node1, node2))

    return new_pairs TODO'''

def create_node_pair_voting(m2m_pairs):
    def add_to_voting(node1, node2):
        if node1 not in node_pair_voting:
            node_pair_voting[node1] = defaultdict(int)

        if node2 not in node_pair_voting:
            node_pair_voting[node2] = defaultdict(int)

        node_pair_voting[node1][node2] += 1
        node_pair_voting[node2][node1] += 1

    node_pair_voting = dict()

    for s1_node, s2_node in m2m_pairs:
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
