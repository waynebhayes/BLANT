#!/pkg/python/3.7.4/bin/python3
import re
import sys
from graph_helpers import *

ORTHO_FILE_PATH = '/home/wayne/src/bionets/SANA/Jurisica/IID/Orthologs.Uniprot.tsv'
SYEAST0_PATH = '../networks/syeast0/syeast0.el'

class SelfOrthos(dict):
    def __setitem__(self, key):
        raise AssertionError

    def __getitem__(self, key):
        return key


def get_avg_node_correctness(all_seeds, g1_to_g2_orthologs):
    if len(all_seeds) == 0:
        return None

    nc_sum = 0

    for gid, align1, align2 in all_seeds:
        num_correct_nodes = 0
        assert len(align1) == len(align2)

        for node1, node2 in zip(align1, align2):
            if is_ortholog(node1, node2, g1_to_g2_orthologs):
                num_correct_nodes += 1

        nc_sum += num_correct_nodes / len(align1)

    return nc_sum / len(all_seeds)

def get_ortho_coverage(all_seeds, g1_to_g2_orthologs):
    orthonodes1 = set()
    orthonodes2 = set()
    g1_orthonodes = set(g1_to_g2_orthologs.keys())
    g2_orthonodes = set(g1_to_g2_orthologs.values())

    for gid, align1, align2 in all_seeds:
        for node in align1:
            if node in g1_orthonodes:
                orthonodes1.add(node)

        for node in align2:
            if node in g2_orthonodes:
                orthonodes2.add(node)

    print(len(orthonodes1), len(orthonodes2))
    return min(len(orthonodes1), len(orthonodes2))

def get_orthoseeds_list(all_seeds_list, s1_to_s2_orthologs, missing_allowed=0):
    orthoseeds_list = []

    for graphlet_id, s1_index, s2_index in all_seeds_list:
        missing_nodes = 0

        assert len(s1_index) == len(s2_index), 's1_index length != s2_index length'

        for m in range(len(s1_index)):
            node1 = s1_index[m]
            node2 = s2_index[m]

            if not is_ortholog(node1, node2, s1_to_s2_orthologs):
                missing_nodes += 1

        if missing_nodes <= missing_allowed:
            orthoseeds_list.append((graphlet_id, s1_index, s2_index))

    return orthoseeds_list

def get_orthopairs_list(node_pairs, s1_to_s2_orthologs):
    orthopairs_list = []

    for node1, node2 in node_pairs:
        if is_ortholog(node1, node2, s1_to_s2_orthologs):
            orthopairs_list.append((node1, node2))

    return orthopairs_list

def is_ortholog(node1, node2, s1_to_s2_orthologs):
    unmarked1 = node1.split('_')[-1]
    unmarked2 = node2.split('_')[-1]

    if type(s1_to_s2_orthologs) is SelfOrthos:
        return unmarked1 == unmarked2
    
    return unmarked1 in s1_to_s2_orthologs and s1_to_s2_orthologs[unmarked1] == unmarked2

def get_g1_to_g2_orthologs(gtag1, gtag2):
    base_gtag1 = gtag1.split('_')[0]
    base_gtag2 = gtag2.split('_')[0]
    g1_is_species = is_species(base_gtag1)
    g2_is_species = is_species(base_gtag2)

    if g1_is_species != g2_is_species:
        raise AssertionError

    if g1_is_species:
        return get_s1_to_s2_orthologs(base_gtag1, base_gtag2)
    else:
        return SelfOrthos()

def get_species_to_index(ortho_file):
    species_to_index = dict()
    species_line = ortho_file.readline().strip()
    species_order = re.split('[\s\t]+', species_line)

    for i, species in enumerate(species_order):
        if species == 'guinea_pig':
            species = 'guineapig'

        species_to_index[species] = i

    return species_to_index

def get_ortholog_nodes(species):
    nodes = list()
    
    with open(ORTHO_FILE_PATH, 'r') as ortho_file:
        species_to_index = get_species_to_index(ortho_file)
        species_pos = species_to_index[species]

        for line in ortho_file:
            line_split = line.strip().split()
            pos_elem = line_split[species_pos]

            if pos_elem == 'guinea_pig':
                pos_elem = 'guineapig'

            if pos_elem == species: # first line
                pass
            else: # other lines
                node = pos_elem

                if node != '0':
                    nodes.append(node)

    return nodes

def get_s1_to_s2_orthologs(species1, species2):
    if 'syeast' in species1 or 'syeast' in species2:
        assert 'syeast' in species1 and 'syeast' in species2, 'syeast must be in both or neither'
        s1_to_s2 = dict()
        nodes = read_in_nodes(SYEAST0_PATH)

        for node in nodes:
            s1_to_s2[node] = node

        return s1_to_s2

    with open(ORTHO_FILE_PATH, 'r') as ortho_file:
        species_to_index = get_species_to_index(ortho_file)
        s1_to_s2 = dict()
        s1_pos = species_to_index[species1]
        s2_pos = species_to_index[species2]

        for line in ortho_file:
            line_split = line.strip().split()

            if line_split[s1_pos] == species1: # first line
                assert line_split[s2_pos] == species2
            else: # other lines
                s1_node = line_split[s1_pos]
                s2_node = line_split[s2_pos]

                if s1_node != '0' and s2_node != '0':
                    s1_to_s2[s1_node] = s2_node

        return s1_to_s2
