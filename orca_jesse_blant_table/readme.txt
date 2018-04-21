orca_jesse_blant_table Readme
This file explains whats in the UpperToLower(k).txt files resulting from running the orca-jesse-blant-table executable
This directory stores both required files and results.
The program generates a table of information for all graphettes from k=3 to 7 and can be used to compare ORCA, Jesse, and BLANT.
This table is sorted in upper row order. 

Columns:
0: CONNECTED:                 1 if the graphette is connected. 0 if disconnected.
1: UPPER_ORDINAL:             The upper ordinal value of the graphette. 
2: UPPER_DECIMAL:             The decimal value of the upper canonical isomorphism of the graphette.
3: LOWER_DECIMAL:             The decimal value of the lower canonical isomorphism of the graphette.
4: LOWER_ORDINAL:             The lower ordinal value of the graphette.
5: NUM_ORBITS:                The number of distinct orbits in the graphette. (From orbit_map[u/l](k).txt)
6: NUM_NODES_FIRST_ORBIT:     The number of nodes of the first orbit type in the graph. (From num_nodes_first_orbit(k).txt)
7: ORCA:                      The ORCA number identifying the graphlet. Prsulj numbering(from GRAAL paper) for k<=5. After, uses lower numbering. Continues form previous k.
8: FIRST_ODV_ORBIT_CON:       The first orbit in a connected orbit degree vector in lower sorting. The sum of previous NUM_ORBITS for connected graphettes. 0 if disconnected
9: FIRST_ODV_ORBIT_ALL:       The first orbit in an orbit degree vector in lower sorting. The sum of previous NUM_ORBITS.
10: FIRST_ODV_ORBIT_CON_FAYE: The first orbit in a connected orbit degree vector in upper sorting. The sum of previous NUM_ORBITS for connected graphettes. 0 if disconnected
11: FIRST_ODV_ORBIT_ALL_FAYE: The first orbit in an orbit degree vector in upper sorting. The sum of previous NUM_ORBITS.

The program fills out columns 0-4 by counting from 0 to the number of canonical nodes for the given k minus one.
This is the UPPER_ORDINAL. It then finds the corresponding UPPER_DECIMAL by looking it up in canon_listu(k).txt.
Then this graphette decimal is converted to lower form, and the canonical is found from that. This is LOWER_ORDINAL.
The LOWER_DECIMAL is found by looking it up in canon_listl(k).txt.

Columns 5 and 6 are loaded from orbit_mapu(k).txt and num_nodes_first_orbit(k).txt respectively.
Column 5 and 0 are used to calculate columns 10 and 11. The algorithm is found in their column description.

The table is then sorted in lower ordering.

Column 5 is then verified from orbit_map(k).txt.
Columns 8 and 9 are found the same way as 10 and 11 but using lower sorting.

If k <=5, column 7 is filled from a manually generated mapping between upper decimal and Prsulj/GRAAL numbering.
If k > 5, the column 7 is filled with the previous total of connected graphettes(graphlets).
If the graphette is disconnected, column 7 is 0.

Finally, the program is resorted in upper ordering.