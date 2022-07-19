#!/pkg/python/3.7.4/bin/python3
from collections import defaultdict

def get_deg_distr(nodes, adj_set):
    deg_distr = defaultdict(int)

    for node in nodes:
        deg = len(adj_set[node])
        deg_distr[deg] += 1

    return deg_distr

def deg_distr_to_str(deg_distr):
    lines = []
    items = list(deg_distr.items())
    items.sort(key=(lambda d: (-d[0], d[1])))
    lines.append('degree\tcount')
    lines.extend([f'{deg}\t{cnt}' for deg, cnt in items])
    return '\n'.join(lines)

def print_deg_distr(deg_distr):
    print(deg_distr_to_str(deg_distr))

def get_seed_nc(seeds, g1_to_g2_ort):
    from ortholog_helpers import is_ortholog

    total_nodes = 0
    weighted_squared_sum = 0

    for gid, nodes1, nodes2 in seeds:
        assert len(nodes1) == len(nodes2), print(gid, nodes1, nodes2, sep='\n')
        seed_size = len(nodes1)
        total_nodes += seed_size
        seed_num_ort = 0
        
        for node1, node2 in zip(nodes1, nodes2):
            if is_ortholog(node1, node2, g1_to_g2_ort):
                seed_num_ort += 1

        weighted_squared_sum += seed_num_ort ** 2 / seed_size # simplified from size * ort ^ 2 / size ^ 2

    weighted_squared_mean = weighted_squared_sum / total_nodes
    return weighted_squared_mean

def get_avg_size(seeds):
    sizes = [len(nodes1) for gid, nodes1, nodes2 in seeds]
    return sum(sizes) / len(sizes)

if __name__ == '__main__':
    gtag1 = 'syeast0'
    gtag2 = 'syeast05'
    g1_to_g2_orts = get_s1_to_s2_orthologs(gtag1, gtag2)
    adj_set1 = read_in_adj_set(get_graph_path(gtag1))
    adj_set2 = read_in_adj_set(get_graph_path(gtag2))
    seeds, _, _ = raw_full_low_param_run(*get_gtag_run_info(gtag1, gtag2))
    orthoseeds = get_orthoseeds_list(seeds, g1_to_g2_orts)
    deg_distr = get_deg_distr(orthoseeds, adj_set1, adj_set2)
    print(len(seeds))
    print_deg_distr(deg_distr)
