#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [OPTIONS] blant.exe M Ks EDs t network.el
PURPOSE: use random samples of k-graphlets from BLANT in attempt to find large communities in network.el.
    blant.exe is the name of the executable BLANT to use (usually just './blant')
    M is the mean number of times each *node* should be touched by a graphlet sample,
    so BLANT will thus take M*(n/k) total samples of k-node graphlets. 
	Ks is a list of graphlet sizes to do BLANT sampling in
	EDs is a list of the edge density thresholds to explore
    t similarity threshold of the output communities
OPTIONS (added BEFORE the blant.exe name)
    -1: exit after printing only one 1 cluster (the biggest one)
    -e: make all the clusters mutually exclusive.
    -o: use only the highest count orbit for neighbors
    -w: cluster count sorted with weights
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

COMMUNITY_MODE=g  # can be g for graphlet (the default if empty), or o for orbit (which uses FAR more memory, like 10-100x)
WEIGHTED=0
ONLY_ONE=0
exclusive=0
SAMPLE_METHOD=-sMCMC #-sNBE
ONLY_BEST_ORBIT=0

while echo "$1" | grep '^-' >/dev/null; do # first argument is an option
    case "$1" in
    -1) ONLY_ONE=1; shift;;
    -o) ONLY_BEST_ORBIT=1; shift;;
    -w) WEIGHTED=1; shift;;
    -e) exclusive=1; die "exclusive not yet supported"; shift;;
    -s) SAMPLE_METHOD="-s$2"; shift 2;;
    -s*) SAMPLE_METHOD="$1"; shift;;
    esac
done

[ $# -lt 6 ] && die "not enough arguments"
[ $# -gt 6 ] && die "too many arguments"

BLANT=$1;
sampleMultiplier=$2;
Ks=(`echo $3 | newlines | sort -nr`); # sort the Ks highest to lowest so the below parallel runs start the higher values of k first
#[ `echo "${Ks[@]}" | wc -w` -eq 1 ] || die "no more multiple K's at the same time"
EDs=($4)
t=$5;
net=$6

numNodes=`awk '{++seen[$1];++seen[$2]}END{print length(seen)}' $net`

[ -x "$BLANT" ] || die "'$BLANT' does not exist or is not an executable"
for k in "${Ks[@]}"; do
    [ "$k" -ge 3 -a "$k" -le 8 ] || die "One k is '$k' but must be between 3 and 8"
done

[ -f "$net" ] || die "network '$net' does not exist"
case "$net" in
*.el) ;;
*) die "network '$net' must be an edgeList file ending in .el";;
esac
DEBUG=false # set to true to store BLANT output

BLANT_EXIT_CODE=0
for k in "${Ks[@]}"; do
    n=`hawk "BEGIN{print int($sampleMultiplier * $numNodes / $k)}"`
    #minEdges=`hawk 'BEGIN{edC='$EDGE_DENSITY_THRESHOLD'*choose('$k',2);rounded_edC=int(edC); if(rounded_edC < edC){rounded_edC++;} print rounded_edC}'`
    # Use MCMC because it gives asymptotically correct concentrations *internally*, and that's what we're using now.
    # DO NOT USE -mi since it will NOT output duplicates, thus messing up the "true" graphlet frequencies/concentrations
    # Possible values: MCMC NBE EBE RES
    [ -f canon_maps/canon_list$k.txt ] || continue
    CMD="$BLANT -k$k -n$n $SAMPLE_METHOD -mc$COMMUNITY_MODE $net" #-e$minEdges
    $CMD > $TMPDIR/blant$k.out &
done

for k in "${Ks[@]}"; do
    wait; (( BLANT_EXIT_CODE += $? ))
done

for edgeDensity in "${EDs[@]}"; do
    for k in "${Ks[@]}"; do
	hawk 'BEGIN{ edC='$edgeDensity'*choose('$k',2); onlyBestOrbit='$ONLY_BEST_ORBIT';
		cMode="'$COMMUNITY_MODE'"; if(cMode=="") cMode=="g"; # graphlet uses FAR less RAM
		ASSERT(cMode=="g" || cMode=="o", "COMMUNITY_MODE must be o or g, not "cMode);
		rounded_edC=int(edC); if(rounded_edC < edC) rounded_edC++;
		minEdges=MAX(rounded_edC, ('$k'-1))
	    }
	    ARGIND==1 && FNR>1 && $2 {canonEdges[FNR-2]=$3}
	    ARGIND==2 && FNR>1 && ((FNR-2) in canonEdges) {for(i=1;i<=NF;i++)orbit2canon[$i]=FNR-2}
	    ARGIND==3 && $3>0 { # ensure the actual count is nonzero
		if(cMode=="o") {
		    orbit=$2; canon=orbit2canon[orbit]; edges=canonEdges[canon]; if(edges<minEdges) next;
		    item=orbit;
		} else if(cMode=="g") {
		    canon=$2; edges=canonEdges[canon]; if(edges<minEdges) next;
		    item=canon;
		}
		if(onlyBestOrbit) item=0;
		Kc[$1][item]+=$3; # increment the cluster count for appropriate item (canon/orbit)
		T[$1]+=$3; # keep total count of *all* items
		# save the neighbors of those cliques that have high edge density for BFS:
		for(j=4;j<=NF;j++) ++neighbors[$1][item][$j];
	    }
	    END {
		for(u in Kc) for(item in Kc[u]) {
		    ORS=" "
		    if('$WEIGHTED') print Kc[u][item]^2/T[u], u, item
		    else            print Kc[u][item], u, item # print the near-clique count and the node
		    for(v in neighbors[u][item]) print v
		    ORS="\n"; print "";
		}
	    }' canon_maps/canon_list$k.txt canon_maps/orbit_map$k.txt $TMPDIR/blant$k.out |
	sort -gr | # > $TMPDIR/cliqs$k.sorted  # sorted near-clique-counts of all the nodes, largest-to-smallest
	hawk 'BEGIN{Srand();OFS="\t"; ID=0;}
	    ARGIND==1{++degree[$1];++degree[$2];edge[$1][$2]=edge[$2][$1]=1} # get the edge list
	    ARGIND==2 && !($2 in count){
		item=$3; count[$2]=$1; node[FNR]=$2; line[$2]=FNR;
		for(i=4; i<=NF; i++) neighbors[$2][$i]=neighbors[$i][$2]=1;
	    }
	    function EdgeCount(v,       edgeHits,u) {
		edgeHits=0;
		for(u in S) if(v in edge[u]){ASSERT(edge[u][v],"edge error"); ++edgeHits;}
		return edgeHits;
	    }
	    function highRelCliqueCount(u, v) { # Heuristic
		if(!(u in count) || count[u]==0) return 1;
		if(!(v in count) || count[v]==0) return 0;
		if (v in count) return count[v]/count[u]>=0.5;
		else return 1/count[u]>=0.5;
	    }
	    function expand(u, origin,    v,oldOrder) {
		if(u in neighbors) {
		    oldOrder=PROCINFO["sorted_in"];
		    PROCINFO["sorted_in"]="randsort";
		    for (v in neighbors[u]) {
			if(!(v in visited) && (!(v in line) || line[v] > line[origin]) && highRelCliqueCount(u, v)) {
			    QueueAdd("Q", v);
			    visited[v]=1;
			}
		    }
		    PROCINFO["sorted_in"]=oldOrder;
		}
	    }
	    END{n=length(degree); # number of nodes in the input network
		cluster[0]=1; delete cluster[0]; # cluster is now explicitly an array, but with zero elements
		numCliques=0;
		QueueAlloc("Q");
		for(start=1; start<=int(FNR*'$edgeDensity'); start++) { # look for a cluster starting on line "start". 
		    if(QueueLength("Q")>0){QueueDelloc("Q");QueueAlloc("Q");} #Ensure the queue is empty
		    origin=node[start];
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
			    if(edgeCount/maxEdges < '$edgeDensity') {
				if(length(S)==Slen){
				    delete S[u]; # no node was removed, so remove this one
				    edgeCount -= newEdgeHits;
				}
				if(++misses > MAX(Slen, n/100)) break; # 1% of number of nodes is a heuristic...
				# keep going until count decreases significantly; duplicate cluster removed in the next awk
				#visited[u]=0;
			    } else
				expand(u, orign);
			} else
			    expand(u, origin);
		    }
		    if(length(S)>'$k') {
			maxEdges=choose(length(S),2);
			++numCliques; printf "%d %d", length(S), edgeCount
			for(u in S) {cluster[numCliques][u]=1; printf " %s", u}
			print ""
			if('$ONLY_ONE') exit;
		    }
		}
	    }' "$net" - | # dash is the output of the above pipe (sorted near-clique-counts)
	sort -nr | # Pablo: Do we need this sort? Eliminating it would allow the below awk to start sooner
	hawk 'BEGIN{ numCliques=0 } # post-process to remove duplicates
	    {
		delete S;
		numNodes=$1
		edgeHits=$2;
		for(i=3;i<=NF;i++) ++S[$i]
		# Pablo: Why is the below line not an error?
		if(length(S)!=numNodes) next;#ASSERT(length(S)==numNodes,"mismatch in numNodes and length(S)");
		add=1;
		for(i=1;i<=numCliques;i++) {
		    same=0;
		    for(u in S) if(u in cluster[i])++same;
		    if(same > length(cluster[i])*'$t'){add=0; break;}
		}
		if(numCliques==0 || add==1) {
		    maxEdges=choose(length(S),2);
		    ++numCliques; 
		    printf "%d %d '$k'",length(S),edgeHits
		    for(u in S) {cluster[numCliques][u]=1; printf " %s", u}
		    print ""
		}
	    }' |
	sort -k 1nr -k 4n > $TMPDIR/subfinal$k$edgeDensity.out & # sort by number of nodes, then by first node in the list
    done
    for k in "${Ks[@]}"; do
	    wait; (( BLANT_EXIT_CODE += $? ))
    done
done

sort -k 1nr -k 4n $TMPDIR/subfinal*.out |
    hawk 'BEGIN{ numCliques=0 } # post-process to remove duplicates
	{
	    delete S; 
	    numNodes=$1
	    edgeHits=$2;
	    k=$3;
	    for(i=4;i<=NF;i++) ++S[$i]
	    if(length(S)!=numNodes) next; #ASSERT(length(S)==numNodes,"mismatch in numNodes and length(S)");
	    add=1;
	    for(i=1;i<=numCliques;i++) {
		same=0;
		for(u in S) if(u in cluster[i])++same;
		if(same > length(cluster[i])*'$t'){ add=0; break;}
	    }
	    if(numCliques==0 || add==1) {
		maxEdges=choose(length(S),2);
		++numCliques; 
		printf "%d nodes, %d of %d edges from k %d (%g%%):", length(S), edgeHits, maxEdges, k, 100*edgeHits/maxEdges
		for(u in S) {cluster[numCliques][u]=1; printf " %s", u}
		print ""
	    }
	}' | sort -k 1nr -k 11n
#set -x
exit $BLANT_EXIT_CODE
