#!/bin/bash

export PATH="`pwd`/libwayne/bin:$PATH" # for hawk
k=5 # only test the 5.el network
G=$k.el # network
export ARGS="-mo -s NBE -n 100000 regression-tests/orcaNumbering/$G" # all ORCA-based BLANT args except -k

# Paste the BLANT outputs for 3,4,5 beside each other (stripping out the leading node names after sorting numerically),
# and then send the output to awk, which will ensure the binary values (zero or nonzero) are identical
paste <(./blant -k 3 $ARGS | sort -n | cut -d' ' -f2-) \
      <(./blant -k 4 $ARGS | sort -n | cut -d' ' -f2-) \
      <(./blant -k 5 $ARGS | sort -n | cut -d' ' -f2-) |
    hawk '{L[ARGIND][FNR]=$0; for(i=1;i<=NF;i++) orbit[ARGIND][FNR][i]=$i}
	END{
	    ASSERT(length(orbit)==2, "wrong number of files");
	    ASSERT(length(orbit[1])==length(orbit[2]), "ARGIND 1 has "length(orbit[1])" lines but 2 has "length(orbit[2]));
	    for(line in orbit[1]){
		ASSERT(length(orbit[1][line])==length(orbit[2][line]),"mismatched line lengths on line "line":"length(orbit[1][line])" vs "length(orbit[2][line]));
		for(col in orbit[1][line])
		    ASSERT(!!orbit[1][line][col] == !!orbit[2][line][col],
			sprintf("mismatched value line %d col %d: %s vs %s:\n%s\n%s",
			    line,col,!!orbit[1][line][col],!!orbit[2][line][col], L[1][line],L[2][line]));
	    }
	}' <(cut -d' ' -f2- regression-tests/orcaNumbering/$k.orca) -
