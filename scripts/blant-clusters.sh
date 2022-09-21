#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
EDGE_DENSITY_THRESHOLD=1.0
USAGE="USAGE: $BASENAME blant.exe k n network.el [ cluster edge density threshold, default $EDGE_DENSITY_THRESHOLD ]
PURPOSE: use n samples of k-graphlets from BLANT in attempt to find large clusters in network.el.  The last argument,
    EDGE_DENSITY_THRESHOLD, is optional and defaults to $EDGE_DENSITY_THRESHOLD."

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
#echo "TMPDIR is $TMPDIR"

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

$BLANT -k$k -n$n -sMCMC -mi "$net" | tee $TMPDIR/blant.out | # produce BLANT index
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
		print Kc[u],u; # print the near-clique count and the node
	    }
	}' <(nl -v -1 canon_maps/canon_list$k.txt) - | # the dash is the BLANT output from -mi run at the top
	    sort -nr | tee $TMPDIR/cliqs.sorted | # sorted near-clique-counts of all the nodes, largest-to-smallest
    hawk 'BEGIN{k='$k';OFS="\t"; ID=0}
	ARGIND==1{edge[$1][$2]=edge[$2][$1]=1} # get the edge list
	ARGIND==2{count[FNR]=$1; node[FNR]=$2}
	function EdgeCount(       edgeCount,u,v) {
	    edgeCount=0;
	    for(u in S){for(v in S) if(u>v && edge[u][v]) ++edgeCount;}
	    return edgeCount;
	}
	END{
	    clique[0]=1; delete clique[0]; # clique is now explicitly an array, but with zero elements
	    numCliques=0;
	    for(start=1; start<=FNR; start++) { # look for a clique starting on line "start"
		delete S; # this will contain the nodes in the current cluster
		lastGood=start;
		S[node[start]]=1;
		for(line=start+1;line<=FNR;line++) {
		    S[node[line]]=1;
		    Slen = length(S);
		    maxEdges=choose(Slen,2);
		    edgeHits = EdgeCount();
		    if(edgeHits/maxEdges < '$EDGE_DENSITY_THRESHOLD') {
			# This is where the greedy algorithm mail fail badly: it is possible that
			# deleting a node currently in S and *keeping* this one may ultimately lead to a larger cluster;
			# we leave this possibility for later implementation. Probably the best possibility is to keep
			# a list of the top X% (X about 80% maybe?) nodes and then use Simulated Annealing to find the
			# biggest clique... but that is much better done in C/C++, not awk.

			numDel=Slen/4 # heuristic
			if(Slen>numDel+3) {
			    # See if removing a recent node or two helps
			    maxEdges1=choose(Slen-1,2);
			    maxEdges2=choose(Slen-2,2);
			    for(del=1; del<numDel; ++del) {
				delete S[node[line-del]]; if(EdgeCount()/maxEdges1 >= '$EDGE_DENSITY_THRESHOLD') break;
				for(del2=del+1; del2<=numDel;++del2) {
				    delete S[node[line-del2]]; if(EdgeCount()/maxEdges2 >= '$EDGE_DENSITY_THRESHOLD') break;
				    ++S[node[line-del2]]; # add it back in
				}
				if(length(S)==Slen-2) break; # removing both nodes helped
				++S[node[line-del]]; # add it back in
			    }
			}
			if(length(S)==Slen) delete S[node[line]]; # no node was removed, so remove this one
			# keep going until count decreases significantly; duplicate cliques removed in the next awk
			if(count[line]/count[lastGood] < 0.5) break;
		    }
		}
		if(length(S)>k) { # now check if it is a subclique of something previously found
		    for(i=1;i<=numCliques;i++) {
			same=0;
			for(u in S) if(u in clique[i]) ++same;
			if(same == length(S)) break;
		    }
		    if(numCliques==0 || same < length(S)) {
			edgeHits=0; maxEdges=choose(length(S),2);
			for(u in S){for(v in S) if(u>v && edge[u][v]) ++edgeHits;}
			++numCliques; printf "%d %d", length(S), edgeHits
			for(u in S) {clique[numCliques][u]=1; printf " %s", u}
			print ""
		    }
		}
	    }
	}' "$net" - | # dash is the output of the above pipe (sorted near-clique-counts)
    sort -nr | # sort the above output by number of nodes in the near-clique
    hawk 'BEGIN{numCliques=0} # post-process to remove duplicates
	{
	    delete S; seenColon=0;
	    numNodes=$1
	    edgeHits=$2;
	    for(i=3;i<=NF;i++) ++S[$i]
	    ASSERT(length(S)==numNodes,"mismatch in numNodes and length(S)");
	    for(i=1;i<=numCliques;i++) {
		same=0;
		for(u in S) if(u in clique[i])++same;
		if(same == length(S)) break;
	    }
	    if(numCliques==0 || same < length(S)) {
		maxEdges=choose(length(S),2);
		++numCliques; printf "%d nodes, %d of %d edges (%g%%):",
		    length(S), edgeHits, maxEdges, 100*edgeHits/maxEdges
		for(u in S) {clique[numCliques][u]=1; printf " %s", u}
		print ""
	    }
	}' | sort -k 1nr -k 8n # sort by number of nodes and then by the first node in the list
