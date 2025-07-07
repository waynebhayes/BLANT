#!/bin/sh
NL='
'
SKIP_EDGES="true" # call with "-e" or -"-a" as only argument to get per-edge details

case "$1" in
-[ae]) SKIP_EDGES=false; shift;;
esac

if "$SKIP_EDGES"; then
    # From https://oeis.org/A000273. Not sure if it includes self-loops.
    echo "k numDirectedCanon"
    echo 1 1 3 16 218 9608 1540944 882033440 1793359192848 13027956824399552 341260431952972580352 32522909385055886111197440 11366745430825400574433894004224 14669085692712929869037096075316220928 70315656615234999521385506555979904091217920 |
	tr ' ' "$NL" | awk '{print NR-1,$0}'
else
    # From https://oeis.org/A052283
    # Triangle read by rows: T(n,k) is the number of unlabeled directed graphs on n nodes with k arcs, k=0..n*(n-1). 
    # Note: includes disconnected graphlets; no self-loops; only up to k=6 nodes
    echo "k	maxEdge	#canon	per_edge_count"
    echo 1 1 1 1 1 1 1 4 4 4 1 1 1 1 5 13 27 38 48 38 27 13 5 1 1 1 1 5 16 61 154 379 707 1155 1490 1670 1490 1155 707 379 154 61 16 5 1 1 1 1 5 17 76 288 1043 3242 8951 21209 43863 78814 124115 171024 207362 220922 207362 171024 124115 78814 43863 21209 8951 3242 1043 288 76 17 5 1 1 |
	hawk '{
	    c=0; # column number from the hugely long echo above (monotonically increases)
	    for(k=0;k<=6;k++) { # graphlet size, 0 through 6 inclusive
		na=MAX(0,k*(k-1)); # na = maximum number of "arcs" (ie., directed edges)=number of columns used for this k
		head=sprintf("%d\t%d",k,na);
		s=0; # sum of canonicals for this k
		sep="\t"; # separator, which changes to a space after first one is printed
		line=""
		for(i=0;i<=na;i++){
		    ++c; # next column
		    s+=$c; # sum
		    line=line sprintf("%s%d",sep,$c);
		    sep=" ";
		}
		printf "%s\t%d%s\n",head,s,line;
	    }
	}'
fi
