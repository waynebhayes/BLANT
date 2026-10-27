#!/bin/sh
# this script takes ORCA output (which we refer to as ODVs, not GDVs, since they are technically _orbit_ count vectors),
# and uses the known overcount factors for each orbit to produce actual _graphlet_ count vectors, or GDVs.
# We do this by picking ONE orbit from each graphlet with a known overcount, and extract that column from ORCA,
# and divide the value by the known overcount.
(
	# First provide the column numbers for the chosen orbits
	echo 1 2 4 5 7 9 10 13 15 16 19 23 25 28 32 35 36 40 44 46 50 52 55 57 60 63 66 69 71 73 | newlines | sed 's/^/orbitCol /'
	# Then provide the overcount factors for the columns selected above
	echo 2 2 3 2 3 4  1  2  4  2  2  4  2  1  2  5  1  1  4  1  3  1  3  1  2  1  1  4  2  5 | newlines | sed 's/^/overcount /'
	# Now spit out whatever file was given on the command line, which must be raw ORCA output for k=5
	cat "$@"
) | # Now pipe those three things--the overcounts, the column selectors, and finally the input--to awk for processing
    awk '
    /^orbitCol/{orbitCol[++numOrbs]=$2;next} # store the column selectors, indexed by "numOrbs"
    /^overcount/{overcount[++oc]=$2;next} # store overcount values, indexed by "oc"

    # Now process the input file(s)
    {
	printf "%d\t", lineNum++;
	space="";
	for(i=1;i<=numOrbs;i++){
	    printf "%s%d", space, $orbitCol[i]/overcount[i];
	    colSum[i]+=$orbitCol[i];
	    space=" ";
	}
	print ""; # end of line
    }

    END{printf "TOTALS\t";
	space="";
	for(i=1;i<=numOrbs;i++) {printf "%s%d", space, colSum[i]/overcount[i]; space=" ";}
	print ""
	}'

# Note both R and awk are indexed by from 1, but the above counts include the edge (g0), so we subtract 1 on output to get ORCA ordering
# Original R code from Sridevi Maharaj:
#counts <- colSums(read.table('/Users/Sridevi/Desktop/RoxannaProject/ORCACounts/orca_con_17_el'))[graphlets]/overcount
