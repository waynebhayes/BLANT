#!/pkg/python/3.7.4/bin/python3
from collections import defaultdict
from general_helpers import *

class Index(defaultdict):
    def __init__(self):
        defaultdict.__init__(self, list)

    def add_entry(self, entry):
        self[entry.get_key()].append(entry)

    def __eq__(self, other):
        if self.keys() != other.keys():
            return False

        for key in self:
            if set(self[key]) != set(other[key]):
                return False

        return True


class IndexEntry:
    def __init__(self, entry_str):
        splitted_str = entry_str.split(' ')
        self._graphlet_id = splitted_str[0] # not an int because we might augment it with bno (base node orbit)
        self._node_arr = splitted_str[1:]
        self._index_of = dict()

        for i, node in enumerate(self._node_arr):
            self._index_of[node] = i

        self._node_set = set(splitted_str[1:])

    def get_key(self):
        return self._graphlet_id

    def get_node_arr(self):
        return self._node_arr

    def node_in(self, node):
        return node in self._node_set

    def index_of(self, node):
        if node in self._index_of:
            return self._index_of[node]
        else:
            return -1

    def get_num_matching(self, other_index):
        num_matching = 0

        for node in self._node_arr:
            if node in other_index._node_set:
                num_matching += 1

        return num_matching

    def get_matching_list(self, other_index):
        matching_list = []

        for i, node in enumerate(self._node_arr):
            if node in other_index._node_set:
                matching_list.append(i)

        return matching_list

    def __str__(self):
        return ' '.join(self._node_arr)

    def __len__(self):
        return len(self._node_arr)

    def __eq__(self, other):
        return set(self._node_arr) == set(other._node_arr)

    def __hash__(self):
        return hash(frozenset(self._node_arr))


# algo None refers to the original algorithm (Henry's)
def get_index_path(gtag, algo='bno', percent=0, orbit=0, alph=True, lDEG=None):
    if lDEG == None:
        if 'syeast' in gtag:
            # lDEG = 3
            lDEG = 2
        else:
            lDEG = 2

    if alph == None: # if alph is none, use old index directory
        index_dir = 'messy_blant_out'
        add_on = ''
    else:
        index_dir = 'blant_out'
        add_on = '-alph' if alph else '-rev'

    if algo == None:
        algo_str = ''
    else:
        algo_str = f'-{algo}'

    return f'{CACHE_BASE_DIR}/{index_dir}/p{percent}-o{orbit}-{gtag}{algo_str}-lDEG{lDEG}{add_on}.out'

def get_notopedge_index_path(gtag, percent=0, orbit=0, lDEG=2):
    return f'{CACHE_BASE_DIR}/special_blant_out/p{percent}-o{orbit}-{gtag}-lDEG{lDEG}-notopedge.out'

def read_in_index(index_path, k):
    with open(index_path, 'r') as index_file:
        index = Index()

        for curr_entry_str in index_file:
            curr_entry_str = curr_entry_str.strip()
            assert len(curr_entry_str.split(' ')) == k + 1, f'the line {curr_entry_str} is not of size k{k}'
            curr_entry = IndexEntry(curr_entry_str)
            index.add_entry(curr_entry)

        return index

def read_in_entry_list(index_path, k):
    with open(index_path, 'r') as index_file:
        entry_list = []

        for curr_entry_str in index_file:
            curr_entry_str = curr_entry_str.strip()
            assert len(curr_entry_str.split(' ')) == k + 1, f'the line {curr_entry_str} is not of size k{k}'
            curr_entry = IndexEntry(curr_entry_str)
            entry_list.append(curr_entry)

        return entry_list

def get_node_distr(index):
    distr = defaultdict(int)

    for entry_list in index.values():
        for entry in entry_list:
            for node in entry.get_node_arr():
                distr[node] += 1

    return distr

if __name__ == '__main__':
    index = read_in_index(get_index_path('syeast0', lDEG=2), 8)
    node_distr = get_node_distr(index)
    print_dict(node_distr)
