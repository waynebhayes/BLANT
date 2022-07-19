#!/bin/python3
from collections import defaultdict

# takes in necessary inputs and settings and returns a list of all found seeds
def find_seeds(s1_index, s2_index, print_progress=True):
    total_pairs_to_process = estimate_total_pairs_to_process(s1_index, s2_index)
    all_seeds_list = []
    percent_printed = 0
    candidate_seeds_processed = 0 # this won't be equal to len(all_seeds_list) because we don't add all candidate seeds

    for graphlet_id, s1_gid_entries in s1_index.items():
        if graphlet_id not in s2_index:
            continue

        s2_gid_entries = s2_index[graphlet_id]

        if len(s1_gid_entries) > 1:
            continue

        if len(s2_gid_entries) > 1:
            continue

        for i in range(0, len(s1_gid_entries)):
            s1_entry_nodes = s1_gid_entries[i].get_node_arr()

            for j in range(0, len(s2_gid_entries)):
                s2_entry_nodes = s2_gid_entries[j].get_node_arr()
                all_seeds_list.append((graphlet_id, tuple(s1_entry_nodes), tuple(s2_entry_nodes)))
                candidate_seeds_processed += 1

                # print
                if print_progress:
                    if candidate_seeds_processed / total_pairs_to_process * 100 > percent_printed:
                        print(f'{percent_printed}% done', file=sys.stderr)
                        percent_printed += 1

    return clean_seeds(all_seeds_list)

def clean_seeds(seeds):
    return list(set(seeds))

def estimate_total_pairs_to_process(s1_index, s2_index):
    total_graphlet_pairs = 0

    for graphlet_id, s1_gid_entries in s1_index.items():
        if graphlet_id in s2_index:
            s2_gid_entries = s2_index[graphlet_id]

            if len(s1_gid_entries) == 1 and len(s2_gid_entries) == 1:
                total_graphlet_pairs += len(s1_gid_entries) * len(s2_gid_entries)

    return total_graphlet_pairs


