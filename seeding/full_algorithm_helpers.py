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
        s1_index = get_patched_index(k, s1_index_path, s1_graph_path)
        s2_index = get_patched_index(k, s2_index_path, s2_graph_path)
        all_seeds_list = find_seeds(10, species1, species2, s1_index, s2_index, get_odv_file_path(species1), get_odv_file_path(species2), SeedingAlgorithmSettings(max_indices=max_indices, sims_threshold=sims_threshold), print_progress=print_progress)
        all_seeds_lists.append(all_seeds_list)

        if print_progress:
            print(f'done with orbit {orbit}')

    return get_combined_seeds_list(all_seeds_lists)

if __name__ == '__main__':
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    species1 = sys.argv[2] if len(sys.argv) > 2 else 'mouse'
    species2 = sys.argv[3] if len(sys.argv) > 3 else 'rat'
    orbits = [int(n) for n in sys.argv[4].split(',')] if len(sys.argv) > 4 else range(2)
    max_indices = int(sys.argv[5]) if len(sys.argv) > 5 else 15
    sims_threshold = float(sys.argv[6]) if len(sys.argv) > 6 else 0.74
    print_progress = sys.argv[7] if len(sys.argv) > 7 else True
    orthoresults, all_results = full_get_patch_pairs_results(k, species1, species2, orbits, max_indices, sims_threshold, print_progress=print_progress)
    print(f'{len(orthoresults)} / {len(all_results)}')
