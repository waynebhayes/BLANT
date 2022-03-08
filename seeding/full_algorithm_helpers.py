#!/bin/python3
# this file runs the entire algorithm from start to finish, hiding all the internal details
from seeding_algorithm_core import *
from node_pair_extraction_helpers import *
from index_helpers import *
from odv_helpers import *
from ortholog_helpers import *
from patch_helpers import *
from validation_helpers import *

def full_get_combined_seeds(k, species1, species2, orbits, max_indices, sims_threshold, print_progress=False):
    all_seeds_lists = []

    for orbit in orbits:
        s1_index_path = get_index_path(species1, orbit=orbit)
        s2_index_path = get_index_path(species2, orbit=orbit)

        if not (validate_index_file(s1_index_path, k) and validate_index_file(s2_index_path, k)):
            print('skipped {species1}-{species2} o{orbit}')
            continue

        all_seeds_list = find_seeds(k, species1, species2, read_in_index(s1_index_path, k), read_in_index(s2_index_path, k), get_odv_file_path(species1), get_odv_file_path(species2), SeedingAlgorithmSettings(max_indices=max_indices, sims_threshold=sims_threshold), print_progress=print_progress)
        all_seeds_lists.append(all_seeds_list)

        if print_progress:
            print(f'done with orbit {orbit}')

    return get_combined_seeds_list(all_seeds_lists)

def full_get_patch_combined_seeds(k, species1, species2, orbits, max_indices, sims_threshold, print_progress=False):
    all_seeds_lists = []

    for orbit in orbits:
        s1_index_path = get_index_path(species1, orbit=orbit)
        s2_index_path = get_index_path(species2, orbit=orbit)

        if not (validate_index_file(s1_index_path, k) and validate_index_file(s2_index_path, k)):
            print('skipped {species1}-{species2} o{orbit}')
            continue

        s1_graph_path = get_graph_path(species1)
        s2_graph_path = get_graph_path(species2)
        all_seeds_list = patch_find_seeds(k, species1, species2, s1_index_path, s2_index_path, s1_graph_path, s2_graph_path, get_odv_file_path(species1), get_odv_file_path(species2), max_indices, sims_threshold, print_progress)
        all_seeds_lists.append(all_seeds_list)

        if print_progress:
            print(f'done with orbit {orbit}')

    return get_combined_seeds_list(all_seeds_lists)

def patch_find_seeds(k, species1, species2, s1_index_path, s2_index_path, s1_graph_path, s2_graph_path, s1_odv_path, s2_odv_path, max_indices, sims_threshold, print_progress=False):
    s1_index = get_patched_index(k, s1_index_path, s1_graph_path)
    s2_index = get_patched_index(k, s2_index_path, s2_graph_path)
    # 10 is the standard size of patched graphlets (patching two 8-node graphlets with 6 nodes in common)
    return find_seeds(10, species1, species2, s1_index, s2_index, s1_odv_path, s2_odv_path, SeedingAlgorithmSettings(max_indices=max_indices, sims_threshold=sims_threshold), print_progress=print_progress)

if __name__ == '__main__':
    k = int(sys.argv[1])
    species1 = sys.argv[2]
    species2 = sys.argv[3]
    s1_index_path = sys.argv[4]
    s2_index_path = sys.argv[5]
    s1_graph_path = sys.argv[6]
    s2_graph_path = sys.argv[7]
    s1_odv_path = sys.argv[8]
    s2_odv_path = sys.argv[9]
    max_indices = int(sys.argv[10])
    sims_threshold = float(sys.argv[11])
    print_progress = sys.argv[12] == 'True'
    all_graphlets = patch_find_seeds(k, species1, species2, s1_index_path, s2_index_path, s1_graph_path, s2_graph_path, s1_odv_path, s2_odv_path, max_indices, sims_threshold, print_progress)
    s1_to_s2_orthologs = get_s1_to_s2_orthologs(species1, species2)
    orthographlets = get_orthographlets_list(all_graphlets, s1_to_s2_orthologs)
    print(f'{len(orthographlets)} / {len(all_graphlets)}')
