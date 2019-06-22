import sys
import numpy as np

num_connected_orbits = {3: 3, 4: 11, 5: 58}
start_index = {3: 1, 4: 4, 5: 15}
exit_code = 0
for k in [3, 4, 5]:
    odv = np.loadtxt("{0}/{1}.odv".format(sys.argv[1], k), delimiter = ' ', usecols=range(1, num_connected_orbits[k]+1))
    orca = np.loadtxt("regression-tests/orcaNumbering/{0}.orca".format(k), delimiter = ' ', usecols=range(start_index[k], start_index[k] + num_connected_orbits[k]))
    if not (np.array_equal(np.array(odv, dtype=bool), np.array(orca, dtype=bool))):
        print("Error: ODV for k={0} doesn't match ORCA".format(k))
        exit_code = 1
sys.exit(exit_code)
