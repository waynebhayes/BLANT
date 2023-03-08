#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
EDGE_DENSITY_THRESHOLD=1.0
USAGE="USAGE: $BASENAME [OPTIONS] blant.exe 'k1 k2...' M network.el [ cluster edge density threshold, default $EDGE_DENSITY_THRESHOLD ] [(k1 k2 ...)]
PURPOSE: use random samples of k-graphlets from BLANT in attempt to find large clusters in network.el.
    blant.exe is the name of the executable BLANT to use (usually just './blant')
    k1 k2...: value(s) of k to use. Multiple values of k can be put in quotes (eg '3 4 5').
    M is the mean number of times each *node* should be touched by a graphlet sample,
	so BLANT will thus take M*(n/k) total samples of k-node graphlets. 
    EDGE_DENSITY_THRESHOLD is optional and defaults to $EDGE_DENSITY_THRESHOLD.
OPTIONS (added BEFORE the blant.exe name)
    -1: exit after printing only one 1 cluster (the biggest one)
    -e: make all the clusters mutually exclusive.
    -sSAMPLE_METHOD: BLANT's sampling method (reasonable choices are MCMC, NBE, EBE, or RES)
    An optional leading integer (with no dash) 'tryHard' [default 0] should not be changed [experimental].
"
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
SAMPLE_METHOD=-sNBE
while echo "$1" | grep '^-' >/dev/null; do # first argument is an option
    case "$1" in
    -1) ONLY_ONE=1; shift;;
    -e) exclusive=1; die "exclusive not yet supported"; shift;;
    -s) SAMPLE_METHOD="$2"; shift 2;;
    -s*) SAMPLE_METHOD="$1"; shift;;
    esac
done

[ $# -lt 4 ] && die "not enough arguments"
#[ $# -gt 6 ] && die "too many arguments"

BLANT=$1;
Ks=(`echo $2 | newlines | sort -nr`); # sort the Ks highest to lowest so the below parallel runs start the higher values of k first
sampleMultiplier=$3 
net=$4; 
if [ $# -eq 5 ]; then
    EDGE_DENSITY_THRESHOLD=$5;
fi

if [ $# -gt 5 ]; then
    EDGE_DENSITY_THRESHOLD=$5; shift 5
fi

numNodes=`newlines < $net | sort -u | wc -l`

[ -x "$BLANT" ] || die "'$BLANT' does not exist or is not an executable"
for k in "${Ks[@]}";
    do
	[ "$k" -ge 3 -a "$k" -le 8 ] || die "One k is '$k' but must be between 3 and 8"
    done

[ -f "$net" ] || die "network '$net' does not exist"
case "$net" in
*.el) ;;
*) die "network '$net' must be an edgeList file ending in .el";;
esac
DEBUG=false # set to true to store BLANT output

for k in "${Ks[@]}";
    do
	n=`hawk 'BEGIN{print int('$sampleMultiplier' * '$numNodes' / '$k')}'`
	edgesCount=`hawk 'BEGIN{edC='$EDGE_DENSITY_THRESHOLD'*choose('$k',2);rounded_edC=int(edC); if(rounded_edC < edC){rounded_edC++;} print rounded_edC}'`
	# DO NOT USE MCMC! Because although MCMC gives asymptotically correct concentrations *internally*, the
	# -mi output will NOT output duplicates, thus messing up the "true" graphlet frequencies/concentrations
	# values: MCMC NBE EBE RES
	echo "[DEBUG=$DEBUG] running: $BLANT -k$k -n$n $SAMPLE_METHOD -mi -e$edgesCount '$net'" >&2
	$BLANT -k$k -n$n $SAMPLE_METHOD -mi -e$edgesCount "$net" > $TMPDIR/blant$k.out & # run them all parallel in the background, outputting to separate files
    done

for k in "${Ks[@]}"; do wait; done

hawk 'BEGIN{Srand(); # fancy random seed
	PROCINFO["sorted_in"]="randsort"; # make for loops go in random order
      }
	{ # run this on ALL input files, not just ARGIND==1
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
			ORS="\n"; print "";
		}
	}' $TMPDIR/blant?.out | # the ? matches all values of k
	sort -nr > $TMPDIR/cliqs.sorted  # sorted near-clique-counts of all the nodes, largest-to-smallest

hawk 'BEGIN{Srand(); # fancy random seed
	PROCINFO["sorted_in"]="randsort"; # make for loops go in random order
	    OFS="\t"; ID=0;
	   }
	ARGIND==1{++degree[$1];++degree[$2];edge[$1][$2]=edge[$2][$1]=1} # get the edge list
	ARGIND==2{count[$2]=$1; node[FNR]=$2; line[$2]=FNR; for(i=3; i<=NF; i++){neighbors[$2][$i]=1;neighbors[$i][$2]=1;}}
	function EdgeCount(v,       edgeHits,u) {
		edgeHits=0;
		for(u in S){if(edge[u][v]) ++edgeHits;}
		return edgeHits;
	}
	function highRelCliqueCount(u, v){ # Heuristic
		if (v in count){
			return count[v]/count[u]>=0.5;
		} else {
			return 1/count[u]>=0.5;
		}
	}

	function expand(u, origin){
		for (v in neighbors[u]){
			if(!(v in visited) && (!(v in line) || line[v] > line[origin]) && highRelCliqueCount(u, v)){
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
		for(start=1; start<=FNR; start++) { # look for a clique starting on line "start". 
			if(QueueLength("Q")>0){QueueDelloc("Q");QueueAlloc("Q");} #Ensure the queue is empty
			origin=node[start];
			if (origin in visited) continue;
			delete S; # this will contain the nodes in the current cluster
			delete visited;
			misses=0; # how many nodes have been skipped because they did not work?
			QueueAdd("Q", origin);
			visited[origin] = 1;
			edgeCount = 0;
			while(QueueLength("Q")>0){
				u = QueueNext("Q");
				newEdgeHits = EdgeCount(u);
				edgeCount += newEdgeHits;
				S[u]=1;
				Slen = length(S);
				if (Slen>1){
					maxEdges = choose(Slen,2);
					if(edgeCount/maxEdges < '$EDGE_DENSITY_THRESHOLD') {
						if(length(S)==Slen){
							delete S[u]; # no node was removed, so remove this one
							edgeCount -= newEdgeHits;
						}
						if(++misses > n/100) break; # 1% of number of nodes is a heuristic...
						# keep going until count decreases significantly; duplicate cliques removed in the next awk
						#visited[u]=0;
					}else{
						expand(u, orign);
					}
				} else {
					expand(u, origin);
				}
			}
			if(length(S)>3) {
				maxEdges=choose(length(S),2);
				++numCliques; printf "%d %d", length(S), edgeCount
				for(u in S) {clique[numCliques][u]=1; printf " %s", u}
				print ""
				if('$ONLY_ONE') exit;
			}
		}
	}' "$net" $TMPDIR/cliqs.sorted  | # dash is the output of the above pipe (sorted near-clique-counts)
    sort -nr | # sort the above output by number of nodes in the near-clique
    hawk 'BEGIN{Srand(); # fancy random seed
	    PROCINFO["sorted_in"]="randsort"; # make for loops go in random order
	    numCliques=0
	} # post-process to remove duplicates
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
