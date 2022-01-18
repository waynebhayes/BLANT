#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
EDGE_DENSITY_THRESHOLD=0.8
USAGE="USAGE: $BASENAME blant.exe k n network.el [ cluster edge density threshold, default $EDGE_DENSITY_THRESHOLD ]
PURPOSE: use n samples of k-graphlets from BLANT in attempt to find large cliques (or more generally clusters) in network.el
    The last argument, EDGE_DENSITY_THRESHOLD, is optional; we stop adding nodes once the edge density under that."

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "$USAGE${NL}FATAL ERROR in $BASENAME:" "$@" >&2; exit 1; }
[ "$BASENAME" == skel ] && die "$0 is a skeleton Bourne Shell script; your scripts should source it, not run it"
echo "$BASENAME" | grep "[ $TAB]" && die "Shell script names really REALLY shouldn't contain spaces or tabs"
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"
warn(){ (echo "WARNING: $@")>&2; }
not(){ if eval "$@"; then return 1; else return 0; fi; }
newlines(){ awk '{for(i=1; i<=NF;i++)print $i}' "$@"; }
parse(){ awk "BEGIN{print $*}" </dev/null; }

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N


#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

[ $# -lt 4 ] && die "not enough arguments"
[ $# -gt 5 ] && die "too many arguments"

BLANT=$1
k=$2
n=$3
net=$4
if [ $# -eq 5 ]; then
    EDGE_DENSITY_THRESHOLD=$5
fi

[ -x "$BLANT" ] || die "'$BLANT' does not exist or is not an executable"
[ "$k" -ge 3 -a "$k" -le 8 ] || die "k is '$k' but must be between 3 and 8"
[ -f "$net" ] || die "network '$net' does not exist"
case "$net" in
*.el) ;;
*) die "network '$net' must be an edgeList file ending in .el";;
esac

$BLANT -k$k -n$n -sMCMC -mi "$net" | # produce BLANT index
    hawk 'BEGIN{k='$k'; want='$EDGE_DENSITY_THRESHOLD'*choose(k,2)} # want = desired minimum number of edges in the k-graphlet
	ARGIND==1{m[$1]=$4} # actual edges in graphlets, from canon_list$k.txt
	ARGIND==2 && m[$1]>=want{
	    for(i=2;i<=NF;i++){
		++Kc[$i]; # increment the near-clique count for each node in the graphlet
		# for(j=i+1;j<=NF;j++){u=MIN($i,$j);v=MAX($i,$j);++Kp[u][v]} # pair count, apparently not needed...
	    }
	}
	END{
	    for(u in Kc){
		print Kc[u],u; # print the near-clique count
	    }
	}' <(nl -v -1 canon_maps/canon_list$k.txt) - | # the dash is the BLANT output from -mi run at the top
	    sort -nr | # now sort the near-clique-counts of all the nodes, largest-to-smallest
    hawk 'BEGIN{k='$k';OFS="\t"; printf "Initial S:"}
	ARGIND==1{edge[$1][$2]=edge[$2][$1]=1} # get the edge list
	ARGIND==2 {
	    S[$2]=1; # take the top nodes according to near-clique count... yes it is that easy.
	    if(FNR<k) printf " %s",$2
	    else if(FNR==k) printf " %s\n",$2
	    else {
		edgeHits=0;
		for(u in S){for(v in S) if(u>v && edge[u][v]) ++edgeHits;}
		maxEdges=choose(FNR,2);
		if(edgeHits/maxEdges < '$EDGE_DENSITY_THRESHOLD') exit;
		printf "adding node %s gives %d nodes with", $2, length(S)
		printf " %d out of a possible %d edges (edge density %.2f%%)\n", edgeHits, maxEdges, 100*edgeHits/maxEdges
	    }
	}' "$net" - # dash is the output of the above pipe (sorted near-clique-counts)
