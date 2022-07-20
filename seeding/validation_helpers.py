#!/bin/python3
from index_helpers import *
from species_helpers import *

def validate_index_file(index_path, k):
    try:
        with open(index_path, 'r') as index_file:
            for entry_str in index_file:
                if len(entry_str.strip().split(' ')) != k + 1:
                    return False
        
            return True
    except FileNotFoundError:
        print(f'file {index_path} not found')
        return False
    except:
        raise

def validate_range(k, species_list, percent_list=[0], orbit_list=[0]):
    for species in species_list:
        for percent in percent_list:
            for orbit in orbit_list:
                is_valid = validate_index_file(get_index_path(species, percent=percent, orbit=orbit), k)

                if not is_valid:
                    print(f'{species} p{percent} o{orbit} is not valid')

if __name__ == '__main__':
    validate_range(8, get_all_species(), [0], list(range(15)))
