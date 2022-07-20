#!/pkg/python/3.7.4/bin/python3
from seeding_algorithm_core import *
from node_pair_extraction_helpers import *
from ortholog_helpers import *
from patch_helpers import *

def full_run(gtag1, gtag2, g1_index_path, g1_graph_path, g2_index_path, g2_graph_path, prox=2, target_num_matching=1):
    k = 8
    s1_to_s2_orthologs = get_s1_to_s2_orthologs(gtag1, gtag2)
    g1_index = get_patched_index(k, g1_index_path, g1_graph_path, prox=prox, target_num_matching=target_num_matching)
    g2_index = get_patched_index(k, g2_index_path, g2_graph_path, prox=prox, target_num_matching=target_num_matching)
    seeds = find_seeds(g1_index, g2_index, print_progress=False)
    g1_to_g2_ort = get_g1_to_g2_orthologs(gtag1, gtag2)
    seed_metrics, extr_metrics = get_all_metrics(seeds, g1_to_g2_ort)
    return (seeds, seed_metrics, extr_metrics)

def get_all_metrics(seeds, g1_to_g2_ort):
    from analysis_helpers import get_avg_size, get_seed_nc

    avg_size = get_avg_size(seeds)
    seed_nc = get_seed_nc(seeds, g1_to_g2_ort)
    seed_metrics = (avg_size, seed_nc)

    all_node_pairs = extract_node_pairs(seeds)
    extr_vol = len(all_node_pairs)
    extr_nc = len(get_orthopairs_list(all_node_pairs, g1_to_g2_ort))
    extr_metrics = (extr_vol, extr_nc)

    return (seed_metrics, extr_metrics)

if __name__ == '__main__':
    gtag1 = sys.argv[1]
    gtag2 = sys.argv[2]
    g1_index_path = sys.argv[3]
    g1_graph_path = sys.argv[4]
    g2_index_path = sys.argv[5]
    g2_graph_path = sys.argv[6]
    seeds, seed_metrics, extr_metrics = full_run(gtag1, gtag2, g1_index_path, g1_graph_path, g2_index_path, g2_graph_path, prox=2, target_num_matching=1)
    print(f'num seeds: {len(seeds)}')
    print(f'avg seed size: {seed_metrics[0]:.3f}')
    print(f'seed node correctness: {seed_metrics[1]:.3f}')
    print(f'num extracted pairs: {extr_metrics[0]}')
    print(f'num correct extracted pairs: {extr_metrics[1]}')
