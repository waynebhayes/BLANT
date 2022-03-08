#!/bin/python3
import sys
from full_algorithm_helpers import *

search_prev_results = dict() # points from threshold to (ortho, all)

def get_result_at_threshold(species1, species2, threshold):
    if threshold in get_result_at_threshold.threshold_results:
        return get_result_at_threshold.threshold_results[threshold]

    k = 8
    orbits = [0]
    max_indices = 10
    result = full_get_combined_seeds(k, species1, species2, orbits, max_indices, threshold)
    get_result_at_threshold.threshold_results[threshold] = result
    return result

# returns highest threshold which isn't a 0
def find_edge(species1, species2):
    # special init
    get_result_at_threshold.threshold_results = dict()

    # vars init
    search_start_step_size = 8
    search_end_step_size = 1
    search_val_start = 100
    search_val_end = 0
    search_step_size = search_start_step_size
    search_val = search_val_start

    # stage 1: go backward
    while True:
        if search_val < search_val_end:
            search_val = search_val_end - 1
            break

        all_pairs = get_result_at_threshold(species1, species2, search_val / 100)

        if len(all_pairs) > 0:
            break

        search_val -= search_step_size

    # stage 2: circle around answer
    if search_val == search_val_start: # special case
        return search_val_start

    # we already know what the first move is, so make it
    search_step_size //= 2
    search_val += search_step_size
    search_step_size //= 2

    while search_step_size >= search_end_step_size:
        all_pairs = get_result_at_threshold(species1, species2, search_val / 100)

        if len(all_pairs) == 0:
            search_val -= search_step_size
        else:
            search_val += search_step_size

        search_step_size //= 2

    all_pairs = get_result_at_threshold(species1, species2, search_val / 100)

    if len(all_pairs) == 0:
        edge = search_val - search_end_step_size
    else:
        edge = search_val

    return edge

def get_final_answer_seeds(species1, species2, print_progress=False):
    edge = find_edge(species1, species2)
    k = 8
    orbits = list(range(15))
    max_indices = 25
    num_under_edge = 10
    best_threshold = (edge - num_under_edge) / 100
    return full_get_seeds_results(k, species1, species2, orbits, max_indices, best_threshold, print_progress)

def get_final_answer_pairs(species1, species2, print_progress=False):
    edge = find_edge(species1, species2)
    k = 8
    orbits = list(range(15))
    max_indices = 25
    max_indices = 3
    num_under_edge = 10
    best_threshold = (edge - num_under_edge) / 100
    return full_get_pairs_results(k, species1, species2, orbits, max_indices, best_threshold, print_progress)

def full_get_patch_pairs_results(k, species1, species2, orbits, max_indices, sims_threshold, print_progress=False):
    combined_seeds = full_get_patch_combined_seeds(k, species1, species2, orbits, max_indices, sims_threshold, print_progress=print_progress)
    node_pairs = extract_node_pairs(combined_seeds)
    s1_to_s2_orthologs = get_s1_to_s2_orthologs(species1, species2)
    orthopairs = get_orthopairs_list(node_pairs, s1_to_s2_orthologs)
    return (orthopairs, node_pairs)

def get_final_results(species1, species2, mode, print_progress=False):
    # settings
    edge = find_edge(species1, species2)
    k = 8
    orbits = list(range(15))
    max_indices = 3
    num_under_edge = 10
    best_threshold = (edge - num_under_edge) / 100

    # run algorithm
    s1_to_s2_orthologs = get_s1_to_s2_orthologs(species1, species2)
    graphlet_pairs = full_get_patch_combined_seeds(k, species1, species2, orbits, max_indices, best_threshold, print_progress=print_progress)

    if mode == 'graphlets':
        ortholog_graphlet_pairs = get_orthographlets_list(graphlet_pairs, s1_to_s2_orthologs)
        return ortholog_graphlet_pairs, graphlet_pairs
    elif mode == 'nodes':
        node_pairs = extract_node_pairs(graphlet_pairs)
        ortholog_node_pairs = get_orthopairs_list(node_pairs, s1_to_s2_orthologs)
        return ortholog_node_pairs, node_pairs

if __name__ == '__main__':
    if get_seeding_dir() == None:
        print('make sure to set the BLANT_SEED_SUPPL envvar')
        quit()

    species1 = sys.argv[1]
    species2 = sys.argv[2]
    mode = sys.argv[3]

    assert mode in {'nodes', 'graphlets'}

    orthoseeds, all_seeds = get_final_results(species1, species2, mode)
    print(f'final results for {species1}-{species2} on mode {mode}: {len(orthoseeds)} / {len(all_seeds)}')
