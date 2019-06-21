import sys
import numpy as np

num_connected_orbits = {3: 3, 4: 11, 5: 58}
start_index = {3: 1, 4: 4, 5: 15}
exit_code = 0
for k in [3, 4, 5]:
    odv = np.loadtxt(f"{sys.argv[1]}/{k}.odv", delimiter = ' ', usecols=range(1, num_connected_orbits[k]+1))
    orca = np.loadtxt(f"regression-tests/orcaNumbering/{k}.orca", delimiter = ' ', usecols=range(start_index[k], start_index[k] + num_connected_orbits[k]))
    if not (np.array_equal(np.array(odv, dtype=bool), np.array(orca, dtype=bool))):
        print(f"Error: ODV for k={k} doesn't match ORCA")
        exit_code = 1
sys.exit(exit_code)
