#!/bin/python3
from collections import defaultdict
from odv_helpers import *
from index_helpers import *
from general_helpers import *

class SeedingAlgorithmSettings:
    def __init__(self, max_indices=15, sims_threshold=0.79, speedup=1):
        self.max_indices = max_indices
        self.sims_threshold = sims_threshold
        self.speedup = speedup

# takes in necessary inputs and settings and returns a list of all found seeds
def find_seeds(k, species1, species2, s1_index=None, s2_index=None, s1_odv_path=None, s2_odv_path=None, settings=SeedingAlgorithmSettings(), print_progress=False):
    if s1_index == None:
        s1_index_path = get_index_path(species1)
        s1_index = read_in_index(s1_index_path, k)

    if s2_index == None:
        s2_index_path = get_index_path(species2)
        s2_index = read_in_index(s2_index_path, k)

    if s1_odv_path == None:
        s1_odv_path = get_odv_file_path(species1)

    if s2_odv_path == None:
        s2_odv_path = get_odv_file_path(species2)
    
    s1_odv_dir = ODVDirectory(s1_odv_path)
    s2_odv_dir = ODVDirectory(s2_odv_path)
    total_pairs_to_process = estimate_total_pairs_to_process(s1_index, s2_index, settings)
    all_seeds_list = []
    percent_printed = 0
    candidate_seeds_processed = 0 # this won't be equal to len(all_seeds_list) because we don't add all candidate seeds

    for graphlet_id, s1_gid_entries in s1_index.items():
        if graphlet_id not in s2_index:
            continue

        s2_gid_entries = s2_index[graphlet_id]

        if len(s1_gid_entries) > settings.max_indices:
            continue

        if len(s2_gid_entries) > settings.max_indices:
            continue

        for i in range(0, len(s1_gid_entries), settings.speedup):
            s1_entry_nodes = s1_gid_entries[i].get_node_arr()

            for j in range(0, len(s2_gid_entries), settings.speedup):
                s2_entry_nodes = s2_gid_entries[j].get_node_arr()
                missing_nodes = 0

                # determine if we want to make it a seed
                if should_be_seed(s1_entry_nodes, s2_entry_nodes, s1_odv_dir, s2_odv_dir, settings.sims_threshold):
                    all_seeds_list.append((graphlet_id, tuple(s1_entry_nodes), tuple(s2_entry_nodes)))

                candidate_seeds_processed += 1

                # print
                if print_progress:
                    if candidate_seeds_processed / total_pairs_to_process * 100 > percent_printed:
                        print(f'{percent_printed}% done', file=sys.stderr)
                        percent_printed += 1
        
    return all_seeds_list

def get_combined_seeds_list(seed_lists):
    final_set = set()

    for seed_list in seed_lists:
        for seed in seed_list:
            final_set.add(seed)

    return list(final_set)

def estimate_total_pairs_to_process(s1_index, s2_index, settings):
    total_graphlet_pairs = 0

    for graphlet_id, s1_gid_entries in s1_index.items():
        if graphlet_id in s2_index:
            s2_gid_entries = s2_index[graphlet_id]

            if len(s1_gid_entries) <= settings.max_indices and len(s2_gid_entries) <= settings.max_indices:
                total_graphlet_pairs += len(s1_gid_entries) * len(s2_gid_entries)

    total_pairs_to_process = int(total_graphlet_pairs / (settings.speedup ** 2))
    return total_pairs_to_process

def should_be_seed(s1_entry_nodes, s2_entry_nodes, s1_odv_dir, s2_odv_dir, threshold):
    sims = []

    for s1_node, s2_node in zip(s1_entry_nodes, s2_entry_nodes):
        s1_odv = s1_odv_dir.get_odv(s1_node)
        s2_odv = s2_odv_dir.get_odv(s2_node)
        sims.append(s1_odv.get_similarity(s2_odv))

    return mean(sims) >= threshold

if __name__ == '__main__':
    seeds = find_seeds(8, 'mouse', 'rat', settings=SeedingAlgorithmSettings(max_indices=3))
    assert_with_prints(len(seeds), 59, 'len(seeds)')
