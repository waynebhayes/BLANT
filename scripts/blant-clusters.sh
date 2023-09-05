#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [OPTIONS] blant_cmd M Ks EDs network.el
PURPOSE: use random samples of k-graphlets from BLANT in attempt to find large communities in network.el.
REQUIRED ARGUMENTS:
    blant_cmd is the raw BLANT command; usually it's just './blant', but can include options (eg './blant -C')
    M is the mean number of times each *node* should be touched by a graphlet sample, so BLANT will thus take
	M*(n/k) total samples of k-node graphlets. 
    Ks is a list of graphlet sizes to do BLANT sampling in
    EDs is a list of the edge density thresholds to explore
OPTIONS (added BEFORE the blant.exe name)
    -D DIR: use existing blant output files in directory DIR
    -1: exit after printing only one 1 cluster (the biggest one)
    -h: use only the highest count graphlet/orbit for neighbors
    -w: cluster count sorted with weights
    -r SEED: use the integer SEED as the random seed
    -m smallest: ignore clusters/cliques/communities with fewer than this number of nodes (default 9)
    -sSAMPLE_METHOD: BLANT's sampling method (reasonable choices are MCMC, NBE, EBE, or RES)
    -o OVERLAP: (default 0.5); the maximum allowed fractional overlap between two clusters;
	clusters that overlap more that this causes either (a) the second one to be discarded, or
	(b) nodes of the second one to be listed as 'nearby' those in the first.
PREDICTION
    -pN: perform edge prediction of up to N edges for each node on the 'periphery' of a communty because it has
	as many as N too few edges into the community to be added; prediction seems to work better with larger
	k, with k=8 working best at high edge densities, though lower k works OK for lower edge densities;
	significantly increasing M also appears to help (eg 10000 or even 100000).
	NOTE: predicted edges are sent to the standard error stream (aka stderr, cerr, Unix file descriptor 2),
	so they do not litter the pipeline. To view them, it's best to run $BASENAME with '2>&1 >/dev/null'
	appended to its command line."

################## SKELETON: DO NOT TOUCH CODE HERE
# check that you really did add a usage message above
USAGE=${USAGE:?"$0 should have a USAGE message before sourcing skel.sh"}
die(){ echo "FATAL ERROR in $BASENAME:" "$@" "${NL}type $BASENAME with no arguments for help" >&2; exit 1; }
usage(){ echo "$USAGE" >&2; exit 1; }
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

[ $# = 0 ] && usage

minClus=9
RANDOM_SEED=
BLANT_FILES="$TMPDIR"
COMMUNITY_MODE=g  # can be g for graphlet (the default if empty), or o for orbit (which uses FAR more memory, like 10-100x)
WEIGHTED=0
ONLY_ONE=0
exclusive=0
SAMPLE_METHOD=-sMCMC #-sNBE
ONLY_BEST_ORBIT=0
OVERLAP=0.5
EDGE_PREDICT=0

# ensure the subfinal sort is always the same... we sort by size since the edge density is (roughly) constant
SUBFINAL_SORT="-k 1nr -k 4n"

while echo "$1" | grep '^-' >/dev/null; do # first argument is an option
    case "$1" in
    -1) ONLY_ONE=1; shift;;
    -h) ONLY_BEST_ORBIT=1; shift;;
    -w) WEIGHTED=1; shift;;
    -s) SAMPLE_METHOD="-s$2"; shift 2;;
    -s*) SAMPLE_METHOD="$1"; shift;;
    -D) BLANT_FILES="$2"; shift 2;;
    -r) RANDOM_SEED="-r $2"; shift 2;;
    -r[0-9]*) RANDOM_SEED="$1"; shift 1;; # allow seed to be same or separate argument
    -m) minClus="$2"; shift 2;;
    -m[0-9]*) minClus="`echo $1|sed 's/^-m//'`"; shift 1;;
    -o) OVERLAP="$2"; shift 2;;
    -o[0-9]*) OVERLAP="`echo $1|sed 's/^-o//'`"; shift 1;;
    -p[0-9]*) EDGE_PREDICT="`echo $1|sed 's/^-p//'`"; shift 1;;
    -*) die "unknown option '$1'";;
    esac
done

[ $# -ne 5 ] && die "expecting exactly 5 mandatory arguments, but you supplied $#"

BLANT_CMD="$1";
BLANT_EXE=`echo $BLANT_CMD | awk '{print $1}'`
sampleMultiplier=$2;
Ks=(`echo $3 | newlines | sort -nr`); # sort the Ks highest to lowest so the below parallel runs start the higher values of k first
#[ `echo "${Ks[@]}" | wc -w` -eq 1 ] || die "no more multiple K's at the same time"
EDs=($4)
net=$5

netCounts=`hawk '{++seen[$1];++seen[$2]}END{n=length(seen);m=NR;printf "%d\t%d\t%g\n", n, m, m/choose(n,2)}' $net`
numNodes=`echo "$netCounts" | cut -f1`
numEdges=`echo "$netCounts" | cut -f2`
meanED=`echo "$netCounts" | cut -f3`

[ -x $BLANT_EXE ] || die "'$BLANT_EXE' does not exist or is not an executable"
for k in "${Ks[@]}"; do
    [ "$k" -ge 3 -a "$k" -le 8 ] || die "One k is '$k' but must be between 3 and 8"
done

[ -f "$net" ] || die "network '$net' does not exist"
case "$net" in
*.el) ;;
*) die "network '$net' must be an edgeList file ending in .el";;
esac
DEBUG=false # set to true to store BLANT output

# Pre-run the BLANTs for each k, and re-use the files for each edge density.
BLANT_EXIT_CODE=0
if [ "$BLANT_FILES" = "$TMPDIR" ]; then
    for k in "${Ks[@]}"; do
	n=`hawk "BEGIN{print int($sampleMultiplier * $numNodes / $k)}"`
	#minEdges=`hawk 'BEGIN{edC='$EDGE_DENSITY_THRESHOLD'*choose('$k',2);rounded_edC=int(edC); if(rounded_edC < edC){rounded_edC++;} print rounded_edC}'`
	# Use MCMC because it gives asymptotically correct concentrations *internally*, and that's what we're using now.
	# DO NOT USE -mi since it will NOT output duplicates, thus messing up the "true" graphlet frequencies/concentrations
	# Possible values: MCMC NBE EBE RES
	[ -f canon_maps/canon_list$k.txt ] || continue
	CMD="$BLANT_CMD $RANDOM_SEED -k$k -n$n $SAMPLE_METHOD -mc$COMMUNITY_MODE $net" #-e$minEdges
	$CMD > $TMPDIR/blant$k.out &
    done
fi

RANDOM_SEED=`echo $RANDOM_SEED | sed 's/^-r *//'` # the "-r" is needed in CMD above but remove it for awk below.

for k in "${Ks[@]}"; do
    wait; (( BLANT_EXIT_CODE += $? ))
done

for edgeDensity in "${EDs[@]}"; do
    for k in "${Ks[@]}"; do
	hawk ' # This first awk runs quite quickly and does not use much RAM
	    BEGIN{ edC='$edgeDensity'*choose('$k',2); onlyBestOrbit='$ONLY_BEST_ORBIT';
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
		# save the neighbors of those clusters that have high edge density for BFS:
		for(j=4;j<=NF;j++) ++neighbors[$1][item][$j];
	    }
	    END {
		for(u in Kc) for(item in Kc[u]) if(item) { # do not use item==0
		    ORS=" "
		    if('$WEIGHTED') print Kc[u][item]^2/T[u], u, item
		    else            print Kc[u][item], u, item # print the cluster membership count, and the node
		    for(v in neighbors[u][item]) print v
		    ORS="\n"; print "";
		}
	    }' canon_maps/canon_list$k.txt canon_maps/orbit_map$k.txt $BLANT_FILES/blant$k.out |
	sort -gr | # > $BLANT_FILES/clus$k.sorted  # sorted cluster-counts of all the nodes, largest-to-smallest
	hawk ' # This second awk can require large amounts of RAM and CPU, eg 10GB and 12-24h for 9wiki-topcats
	    BEGIN{if("'$RANDOM_SEED'") srand('$RANDOM_SEED');else Srand();OFS="\t"; ID=0; edgePredict='"$EDGE_PREDICT"';}
	    ARGIND==1{++degree[$1];++degree[$2];edge[$1][$2]=edge[$2][$1]=1} # get the edge list
	    ARGIND==2 && !($2 in count){ # are we really SURE we should take only the first occurence?
		item=$3; count[$2]=$1; node[FNR]=$2; line[$2]=FNR;
		for(i=4; i<=NF; i++) neighbors[$2][$i]=neighbors[$i][$2]=1;
	    }

	    function EdgesIntoS(v,       edgeHits,u) { # count edges from v into S, not including v (which can be is S)
		edgeHits=0;
		for(u in S) if(u!=v && (v in edge[u])){
		    ASSERT(edge[v][u],"edge error");
		    ++edgeHits;
		}
		return edgeHits;
	    }
	    function highRelClusCount(u, v) { # Heuristic
		if(!(u in count) || count[u]==0) return 1;
		if(!(v in count) || count[v]==0) return 0;
		if (v in count) return count[v]/count[u]>=0.5;
		else return 1/count[u]>=0.5;
	    }
	    function AppendNeighbors(u, origin,    v,oldOrder) {
		oldOrder=PROCINFO["sorted_in"];
		PROCINFO["sorted_in"]="randsort";
		if(u in neighbors) {
		    for (v in neighbors[u]) {
			if(!(v in visitedQ) && (!(v in line) || line[v] > line[origin]) && highRelClusCount(u, v)) {
			    QueueAdd("Q", v);
			    visitedQ[v]=1;
			}
		    }
		}
		PROCINFO["sorted_in"]=oldOrder;
	    }
	    END{n=length(degree); # number of nodes in the input network
		# make cluster and started empty sets; started is a list of nodes we should NOT start a new BFS on
		delete cluster; cluster[0]=1; delete cluster[0];
		delete started; started[0]=1; delete started[0];
		numClus=0;
		QueueAlloc("Q");
		for(start=1; start<=int(FNR*'$edgeDensity'); start++) { # look for a cluster starting on line "start". 
		    origin=node[start];
		    if(started[origin]) continue;
		    started[origin]=1;
		    MakeEmptySet(S);
		    MakeEmptySet(visitedQ);
		    misses=0; # how many nodes have been skipped because they did not work?
		    if(QueueLength("Q")>0){QueueDelloc("Q");QueueAlloc("Q");} #Ensure the queue is empty
		    QueueAdd("Q", origin);
		    visitedQ[origin] = 1;
		    edgeCount = 0;
		    while(QueueLength("Q")>0){
			u = QueueNext("Q");
			ASSERT(!(u in S),"u in S error");
			newEdgeHits = EdgesIntoS(u);
			edgeCount += newEdgeHits;
			S[u]=1;
			Slen = length(S);
			# CAREFUL: calling InducedEdges(S) every QueueNext() is VERY expensive; uncommenting the line below
			# slows the program by more than 100x (NOT an exaggeration)
			# WARN(edgeCount == InducedEdges(edge,S),"Slen "Slen" edgeCount "edgeCount" Induced(S) "InducedEdges(edge,S));
			Sorder[Slen]=u; # no need to delete this element if u fails because Slen will go down by 1
			if (Slen>1){
			    maxEdges = choose(Slen,2);
			    if(edgeCount/maxEdges < '$edgeDensity') {
				delete S[u]; # u drops the edge density too low, so nuke it
				if(edgePredict) {
				    if((edgeCount+edgePredict)/maxEdges >= '$edgeDensity') {
					# just ONE edge would put us over the threshold; predict it should exist
					numPredict=0;
					for(uu in S)if(!(u in edge[uu]))
					    printf "predictEdge-%d\t%s\t%s\n",++numPredict, u,uu > "/dev/stderr"
					WARN(numPredict<=edgePredict,"should only get "edgePredict" new edge predictions but got "numPredict);
				    }
				}
				edgeCount -= newEdgeHits;
				if(++misses > MAX(Slen, n/100)) break; # 1% of number of nodes is a heuristic...
				# keep going until count decreases significantly; remove duplicate clusters in next awk
			    } else {
				misses=0;
				AppendNeighbors(u, orign);
			    }
			} else
			    AppendNeighbors(u, origin);
		    }
		    if(QueueLength("Q")==0) {
			# post-process to remove nodes that have too low degree compared to the norm. We use two criteria:
			# 1) the in-cluster degree is more than 3 sigma below the mean, or
			# 2) the in-cluster degree is less than 1/3 the mode of the in-cluster degree distribution.
			# The latter was added in response to our performance on the LFR graphs, but it does not appear
			# to hurt performance anywhere else.
			StatReset(""); delete degFreq;
			tmpEdge = InducedEdges(edge,S, degreeInS);
			ASSERT(tmpEdge == edgeCount, "edgeCount "edgeCount" disagrees(1) with InducedEdges of "tmpEdge);
			for(v in S) { StatAddSample("", degreeInS[v]); ++degFreq[degreeInS[v]];}
			maxFreq=degMode=0;
			for(d=StatMax("");d>=StatMin("");d--)
			    if(d in degFreq && degFreq[d] >= maxFreq) { # use >= to extract SMALLEST mode
				maxFreq=degFreq[d]; degMode=d
			}
			#printf "maxFreq %d mode %d\n", maxFreq, degMode > "/dev/stderr"
			#printf "start %d |S|=%d mean %g stdDev %g:", start, _statN[""], StatMean(""), StatStdDev("") > "/dev/stderr"
			for(v in S) if(degreeInS[v] < StatMean("") - 3*StatStdDev("") || degreeInS[v] < degMode/3) {
			    #printf " %s(%d)", v, degreeInS[v] > "/dev/stderr";
			    delete S[v];
			    edgeCount -= EdgesIntoS(v);
			}
		    }
		    Slen=length(S);
		    #printf " final |S|=%d\n",Slen > "/dev/stderr";
		    if(Slen>='$minClus') {
			maxEdges=choose(Slen,2);
			++numClus; printf "%d %d", Slen, edgeCount
			StatReset("");
			tmpEdge=InducedEdges(edge,S, degreeInS);
			ASSERT(tmpEdge == edgeCount, "edgeCount "edgeCount" disagrees(2) with InducedEdges of "tmpEdge);
			for(u in S) {cluster[numClus][u]=1; printf " %s", u; StatAddSample("", degreeInS[u]);}
			#printf "final |S|=%d mean %g stdDev %g min %d max %d:\n", _statN[""], StatMean(""), StatStdDev(""), StatMin(""), StatMax("") > "/dev/stderr"
			print "";
			if('$ONLY_ONE') exit;
			# now mark as "started" the first half of S that was filled... or more generally, some fraction.
			# We choose the top (1-OVERLAP) fraction because OVERLAP is meant to be more stringent as it approaches
			# 1, which in the case of THIS loop means that if OVERLAP is close to 1, we want to eliminate a
			# SMALLER proportion of future BFS starts, so if OVERLAP=1 we eliminate NOTHING in this loop.
			for(i=1;i<Slen*(1-'$OVERLAP');i++) started[Sorder[i]]=1;
		    }
		}
	    }' "$net" - | # dash is the output of the above pipe (sorted cluster counts)
	sort -nr |
	hawk ' # This third awk attempts to remove duplicate clusters... not sure how intensive it is in RAM and CPU (yet)
	    BEGIN{ numClus=0 } # post-process to only EXACT duplicates (more general removal later)
	    {
		delete S;
		numNodes=$1
		edgeHits=$2;
		for(i=3;i<=NF;i++) ++S[$i]
		WARN(length(S)==numNodes,"mismatch(1) in numNodes and length(S)");
		add=1;
		for(i=1;i<=numClus;i++) {
		    same=0;
		    for(u in S) if(u in cluster[i])++same;
		    if(same == length(cluster[i])){add=0; break;} # only eliminate EXACT duplicates at this stage
		}
		if(add) {
		    ++numClus; edges[numClus]=edgeHits;
		    for(u in S) ++cluster[numClus][u]
		}
	    }
	    END{
		for(i=1;i<=numClus;i++) {
		    maxEdges=choose(length(cluster[i]),2);
		    printf "%d %d '$k'",length(cluster[i]),edges[i]
		    for(u in cluster[i]) printf " %s", u
		    print ""
		}
	    }
	    ' | # sort by number of nodes, then by first node in the list:
	sort $SUBFINAL_SORT > $TMPDIR/subfinal$k$edgeDensity.out & # sort by number of nodes, then by first node in the list
    done
    for k in "${Ks[@]}"; do
	    wait; (( BLANT_EXIT_CODE += $? ))
    done
done

# we can use the --merge option because the subfinal files are already sorted in SUBFINAL_SORT order
sort --merge $SUBFINAL_SORT $TMPDIR/subfinal*.out |
    hawk 'BEGIN{ numClus=0 } # post-process to remove/merge duplicates
	{
	    delete S; 
	    numNodes=$1
	    edgeHits=$2;
	    k=$3;
	    for(i=4;i<=NF;i++) ++S[$i]
	    WARN(length(S)==numNodes,"mismatch(2) in numNodes and length(S)");
	    add=1;
	    for(i=1;i<=numClus;i++) {
		same=0;
		for(u in S) if(u in cluster[i])++same;
		if(same > MAX(length(S), length(cluster[i]))*'$OVERLAP'){ add=0; break;}
	    }
	    # Skipping an entire cluster with overlap sucks, but merging does not work either. Here is a better way:
	    # The idea is to provide ALTERNATES.  So for example if 2 clusters of size k overlap in all but one node,
	    # and have the same total number of edges, then the 2 nodes that are OUTSIDE the overlap form "alternate"
	    # forms of an identical cluster, at least in quality. (It is likely that they effectively form identical
	    # graphlets with a large k, and the 2 nodes occupy the same orbit in the graphlet). So if node A can be
	    # swapped out with node B to get an identical quality cluster, then that could be listed as "one" node
	    # A#B in the output, meaning A and B can we swapped and you still get the same number of nodes and edges,
	    # while adding both does not work. Doing this correctly requires some thought, so it will have to wait.

	    if(add) {
		++numClus; edges[numClus]=edgeHits; kk[numClus]=k;
		for(u in S) ++cluster[numClus][u];
	    }
	}
	END {
	    PROCINFO["sorted_in"]="@ind_num_asc"; # print nodes in numerical ascending order
	    #PROCINFO["sorted_in"]="@ind_str_asc"; # print nodes in string ascending order
	    for(i=1;i<=numClus;i++) {
		maxEdges=choose(length(cluster[i]),2);
		printf "%d nodes, %d of %d edges from k %d (%g%%):", length(cluster[i]), edges[i], maxEdges, kk[i], 100*edges[i]/maxEdges
		for(u in cluster[i]) printf " %s", u
		print ""
	    }
	}' | sort -k 1nr -k 11n
#set -x
exit $BLANT_EXIT_CODE
