#!/bin/sh

hawk '{
	if     (NF==15 || NF==16){ k=4; nOrbits=15; nGraphlets= 9; hasNode=(NF>15);}
	else if(NF==73 || NF==74){ k=5; nOrbits=73; nGraphlets=30; hasNode=(NF>73);}
	else ASSERT(false, "wrong number of columns");
	for(i=0;i<nOrbits;i++)
	    orbit_count[i] += $(hasNode + i+1) # +1 since awk numbers columns from 1 but Python and ORCA from 0
    }
    END{
	graphlet_count[0] = orbit_count[0]/2

	graphlet_count[1] = orbit_count[2]
	graphlet_count[2] = orbit_count[3]/3

	graphlet_count[3] = orbit_count[4]/2
	graphlet_count[4] = orbit_count[7]
	graphlet_count[5] = orbit_count[8]/4
	graphlet_count[6] = orbit_count[9]
	graphlet_count[7] = orbit_count[12]/2
	graphlet_count[8] = orbit_count[14]/4

	if(k==5) {
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
	    graphlet_count[21] = orbit_count[52]
	    graphlet_count[22] = orbit_count[55]/2
	    graphlet_count[23] = orbit_count[58]
	    graphlet_count[24] = orbit_count[61]
	    graphlet_count[25] = orbit_count[62]
	    graphlet_count[26] = orbit_count[65]
	    graphlet_count[27] = orbit_count[69]
	    graphlet_count[28] = orbit_count[70]/2
	    graphlet_count[29] = orbit_count[72]/5
	}
	printf "ORCA numbering:"
	for(i=0;i<nGraphlets;i++) printf(" %d",graphlet_count[i]); print "";

	canon2gID[3][3]=1;
	canon2gID[3][7]=2;

	canon2gID[4][7]=4;
	canon2gID[4][13]=3;
	canon2gID[4][15]=6;
	canon2gID[4][30]=5;
	canon2gID[4][31]=7;
	canon2gID[4][63]=8;

	if(k==5) {
	    canon2gID[5][15]=11;
	    canon2gID[5][29]=10;
	    canon2gID[5][31]=14;
	    canon2gID[5][58]=9;
	    canon2gID[5][59]=12;
	    canon2gID[5][62]=16;
	    canon2gID[5][63]=17;
	    canon2gID[5][126]=20;
	    canon2gID[5][127]=22;
	    canon2gID[5][185]=13;
	    canon2gID[5][187]=19;
	    canon2gID[5][191]=23;
	    canon2gID[5][207]=18;
	    canon2gID[5][220]=15;
	    canon2gID[5][221]=21;
	    canon2gID[5][223]=24;
	    canon2gID[5][254]=25;
	    canon2gID[5][255]=26;
	    canon2gID[5][495]=27;
	    canon2gID[5][511]=28;
	    canon2gID[5][1023]=29;
	}

	print "\nBLANT numbering:"
	PROCINFO["sorted_in"]="@ind_num_asc"; # iterate through arrays numerically by index
	for(k in canon2gID) {
	    printf "k = %d", k;
	    for(canon in canon2gID[k]) printf " %d", graphlet_count[canon2gID[k][canon]]
	    print "";
	}
    }' "$@"
