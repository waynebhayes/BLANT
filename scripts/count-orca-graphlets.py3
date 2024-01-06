import sys

tot = 29

orbit_count = []

with open("/tmp/FUCKING-PYTHON-INPUT.SHIT" , "r") as f:
    lines = f.readlines()

    ans = [0 for num in lines[0].split()]

    for line in lines:
        ans = [ans[j] + int(num) for j, num in enumerate(line.split())]

    orbit_count = [int(num) for num in ans[1:]]

graphlet_count = [0 for i in range(tot+1)]
concentration = [1 for i in range(tot+1)]


graphlet_count[0] = orbit_count[0]/2

graphlet_count[1] = orbit_count[2]
graphlet_count[2] = orbit_count[3]/3

graphlet_count[3] = orbit_count[4]/2
graphlet_count[4] = orbit_count[7]
graphlet_count[5] = orbit_count[8]/4
graphlet_count[6] = orbit_count[9]
graphlet_count[7] = orbit_count[12]/2
graphlet_count[8] = orbit_count[14]/4

graphlet_count[9] = orbit_count[17]
graphlet_count[10] = orbit_count[21]
graphlet_count[11] = orbit_count[23]
graphlet_count[12] = orbit_count[25]
graphlet_count[13] = orbit_count[30]
graphlet_count[14] = orbit_count[33]
graphlet_count[15] = orbit_count[34]/5
graphlet_count[16] = orbit_count[38]
graphlet_count[17] = orbit_count[39]
graphlet_count[18] = orbit_count[44]
graphlet_count[19] = orbit_count[47]
graphlet_count[20] = orbit_count[50]/2
graphlet_count[21] = orbit_count[51]
graphlet_count[22] = orbit_count[55]/2
graphlet_count[23] = orbit_count[58]
graphlet_count[24] = orbit_count[61]
graphlet_count[25] = orbit_count[62]
graphlet_count[26] = orbit_count[65]
graphlet_count[27] = orbit_count[69]
graphlet_count[28] = orbit_count[70]/2
graphlet_count[29] = orbit_count[72]/5

print("ORCA numbering:", graphlet_count)

canon2gID = {
    3: {
        3: 1,
        7: 2
    },
    4: {
        7: 4,
        13: 3,
        15: 6,
        30: 5,
        31: 7,
        63: 8
    },
    5: {
        15: 11,
        29: 10,
        31: 14,
        58: 9,
        59: 12,
        62: 16,
        63: 17,
        126: 20,
        127: 22,
        185: 13,
        187: 19,
        191: 23,
        207: 18,
        220: 15,
        221: 21,
        223: 24,
        254: 25,
        255: 26,
        495: 27,
        511: 28,
        1023: 29
    }
}

print("\nBLANT numbering:")
for k in canon2gID:
    print("k =",k, [ graphlet_count[canon2gID[k][canon]] for canon in canon2gID[k]])

