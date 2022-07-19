import sys
from statistics import mean
from graph_helpers import *
from general_helpers import *
from species_helpers import *

def get_odv_file_path(species):
    return f'{get_seeding_dir()}/odv/{species_to_full_name(species)}.orca4'

class ODVDirectory:
    # file format: every line has node name, followed by orbit counts, separated by spaces
    # NODENAME 23 1 250 37 4 0 ...
    def __init__(self, fname):
        self._directory = dict()

        for line in open(fname, 'r'):
            line_split = line.strip().split()
            node = line_split[0]
            odv_list = [int(s) for s in line_split[1:]]
            odv = ODV(odv_list)
            self._directory[node] = odv

    def get_odv(self, node):
        return self._directory[node]

    def get_nodes(self):
        return self._directory.keys()

    def __str__(self):
        return '\n'.join([f'{node}: {odv}' for node, odv in self._directory.items()])


class ODV:
    def __init__(self, odv_list):
        self._odv_list = odv_list

    def get_similarity(self, other):
        return mean([self._get_single_orbit_similarity(m1, m2) for m1, m2 in zip(self._odv_list, other._odv_list)])

    def get_odv_val(self, num):
        return self._odv_list[num]

    def __str__(self):
        return ' '.join([str(n) for n in self._odv_list])

    @staticmethod
    def _get_single_orbit_similarity(m1, m2):
        return 1 if m1 == m2 == 0 else min(m1, m2) / max(m1, m2)
