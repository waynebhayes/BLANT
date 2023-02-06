#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
EDGE_DENSITY_THRESHOLD=1.0
USAGE="USAGE: $BASENAME [tryHard] [-1] [-e] blant.exe k M network.el [ cluster edge density threshold, default $EDGE_DENSITY_THRESHOLD ]
PURPOSE: use random samples of k-graphlets from BLANT in attempt to find large clusters in network.el.
    M is the mean number of times each *node* should be touched by a graphlet sample; thus, BLANT will take M*(n/k) total
    samples of k-node graphlets; we recommend experimenting with a minimum of M=100, and going as high as necessary (over
    1000 is not uncommon) for reliable results.
    The last argument, EDGE_DENSITY_THRESHOLD, is optional and defaults to $EDGE_DENSITY_THRESHOLD.
    An optional leading integer (with no dash) 'tryHard' [default 0] should not be changed [experimental].
    The option '-1' means 'exit after printing only one 1 cluster--the top one'.
    The option '-e' means make all the clusters mutually exclusive."

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
TMPDIR=`mktemp -d ${LOCAL_TMP:-"/tmp"}/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
#echo "TMPDIR is $TMPDIR"

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

tryHard=0
case "$1" in
[0-9]*) tryHard=$1; shift;;
esac

ONLY_ONE=0
exclusive=0
while echo "$1" | grep '^-' >/dev/null; do # first argument is an option
    case "$1" in
    -1) ONLY_ONE=1; shift;;
    -e) exclusive=1; die "exclusive not yet supported"; shift;;
    esac
done

[ $# -lt 4 ] && die "not enough arguments"
[ $# -gt 5 ] && die "too many arguments"

BLANT=$1
k=$2
sampleMultiplier=$3
net=$4
if [ $# -eq 5 ]; then
    EDGE_DENSITY_THRESHOLD=$5
fi

numNodes=`newlines < $net | sort -u | wc -l`
n=`expr $sampleMultiplier '*' $numNodes / $k`

[ -x "$BLANT" ] || die "'$BLANT' does not exist or is not an executable"
[ "$k" -ge 3 -a "$k" -le 8 ] || die "k is '$k' but must be between 3 and 8"
[ -f "$net" ] || die "network '$net' does not exist"
case "$net" in
*.el) ;;
*) die "network '$net' must be an edgeList file ending in .el";;
esac
DEBUG=false # set to true to store BLANT output
edgesCount=`hawk 'BEGIN{edC='$EDGE_DENSITY_THRESHOLD'*choose('$k',2);rounded_edC=int(edC); if(rounded_edC < edC){rounded_edC++;} print rounded_edC}'`
echo "[DEBUG=$DEBUG] running: $BLANT -k$k -n$n -sMCMC -mi '$net' -e$edgesCount " >&2

$BLANT -k$k -n$n -sMCMC -mi "$net" -e$edgesCount  | if $DEBUG; then tee $TMPDIR/blant.out; else cat; fi |
    hawk 'BEGIN{k='$k'} 
		ARGIND==1{
			for(i=2;i<=NF;i++){
				++Kc[$i]; # increment the near-clique count for each node in the graphlet
				for(j=2;j<=NF;j++){ # saving the neighbors of those cliques that have high edge density for BFS
					if (j==i) continue;
					neighbors[$i][$j] = 1
				}
			}
		}
		END{
			for(u in Kc){
				ORS=" "
				print Kc[u], u # print the near-clique count and the node
				for (v in neighbors[u]){
					print v
				}
				ORS="\n"; print;
			}
		}' - | # the dash is the BLANT output from -mi run at the top
	    sort -nr | tee $TMPDIR/cliqs.sorted | # sorted near-clique-counts of all the nodes, largest-to-smallest
    hawk 'BEGIN{k='$k';OFS="\t"; ID=0}
		ARGIND==1{++degree[$1];++degree[$2];edge[$1][$2]=edge[$2][$1]=1} # get the edge list
		ARGIND==2{count[$2]=$1; node[FNR]=$2; line[$2]=FNR; for(i=3; i<=NF; i++){neighbors[$2][$i]=1;}}
		function EdgeCount(       edgeCount,u,v) {
			edgeCount=0;
			for(u in S){for(v in S) if(u>v && edge[u][v]) ++edgeCount;}
			return edgeCount;
		}
		function expand(u){
			for (v in neighbors[u]){
				if((v in line) && !(v in visited)){
					QueueAdd("Q", v);
					visited[v]=1;
				}
			}
			return;
		}
		END{n=length(degree); # number of nodes in the input network
			clique[0]=1; delete clique[0]; # clique is now explicitly an array, but with zero elements
			numCliques=0;
			QueueAlloc("Q");
			for(start=1; start<=FNR; start++) { # look for a clique starting on line "start"
				while(QueueLength("Q")>0){QueueNext("Q");} #Ensure the queue is empty
				delete S; # this will contain the nodes in the current cluster
				delete visited;
				origin=node[start];
				misses=0; # how many nodes have been skipped because they did not work?
				QueueAdd("Q", origin);
				visited[origin] = 1;
				while(QueueLength("Q")>0){
					u = QueueNext("Q");
					S[u]=1;
					Slen = length(S);
					if (Slen>1){
						maxEdges = choose(Slen,2);
						edgeHits = EdgeCount();
						if(edgeHits/maxEdges < '$EDGE_DENSITY_THRESHOLD') {
							if(++misses > n/100) break; # 1% of number of nodes is a heuristic...
							if(length(S)==Slen) delete S[u]; # no node was removed, so remove this one
							# keep going until count decreases significantly; duplicate cliques removed in the next awk
							if(count[u]/count[origin] < 0.5) break;
						}else{
							expand(u);
						}
					} else {
						expand(u);
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
						if('$ONLY_ONE') exit;
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
