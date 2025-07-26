#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
export BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
export LINENO
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [OPTIONS] blant_cmd Ks EDs network.el
PURPOSE: use random samples of k-graphlets from BLANT in attempt to find large communities in network.el.
REQUIRED ARGUMENTS:
    blant_cmd is the raw BLANT command; usually it's just './blant', but can include options (eg './blant -C')
    Ks is a list of graphlet sizes to do BLANT sampling in
    EDs is a list of the edge density thresholds to explore
OPTIONS (added BEFORE the blant.exe name)
    -D DIR: use existing blant output files in directory DIR
    -T DIR: use temp space here rather than /tmp
    -A: Use Sweta Jain's Turan's Shadow code to estimate absolute counts
    -C: run on the complement graph, G' (careful, this can cause long run times)
    -b: exit after printing only the biggest cluster
    -1: only take the first occurence of a node and its neighbors (faster but may get worse result)
    -h: use only the highest count graphlet/orbit for neighbors
    -S: cluster count sorted by the normalized square
    -w: networks are edge-weighted; pass -w to BLANT to use weighted graphlets
    -v: verbose output
    -nm: 'no members': suppress list of nodes in the community
    -pD: tell BLANT to estimate graphlet counts with D digits of precision (default 1)
    -nN: tell BLANT to take exactly N samples (mutually exclusive with -p above)
    -r SEED: use the integer SEED as the random seed
    -emt value: minimum edge mean threshold: the smallest average-over-edgeWeights for a cluster to be included [default 0.5]
    -m smallest: if m>=1, ignore clusters/cliques/communities with fewer than this number of nodes; otherwise if m<1, treat
       it as a p-value and set the smallest cluster to whatever size gives the number of edges a p-value less than this;
       default behaviour is to demand a cluster whose number of edges has p-value < 1/(n choose 2)^2 (heuristic).
    -M largest: stop building a community if it reaches this size (default: no limit)
    -sSAMPLE_METHOD: BLANT's sampling method (reasonable choices are MCMC, NBE, EBE, or RES)
    -o OVERLAP: (default 0.5); the maximum allowed fractional overlap between a new cluster and any one previously discovered
	cluster (though it can overlap by more than this with MULTIPLE previous clusters). If it overlaps more than this with
	any ONE previous cluster, it is discarded.
PREDICTION
    -PN: perform edge prediction of up to N edges for each node on the 'periphery' of a communty because it has
	as many as N too few edges into the community to be added; prediction seems to work better with larger
	k, with k=8 working best at high edge densities, though lower k works OK for lower edge densities.
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

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

[ $# = 0 ] && usage

minClusArg=0
maxClus=2147483647 # 2^31-1
minEdgeMean=0.5
RANDOM_SEED=
COMMUNITY_MODE=g  # can be g for graphlet (the default if empty), or o for orbit (which uses FAR more memory, like 10-100x)
SQR_NORM=0
WEIGHTED=''
W=''
ONLY_BIGGEST=0
exclusive=0
SAMPLE_METHOD=-sEBE #MCMC #-sNBE
ONLY_BEST_ORBIT=0
ONLY_FIRST=0
OVERLAP=0.5
EDGE_PREDICT=0
COMPLEMENT=''
TURAN=false
DEBUG=false # set to true to store BLANT output
VERBOSE=0
QUIET='' #-qq
PRECISION=-p1L
PRINT_MEMBERS=1
DENSITY_LEEWAY=0.95 # factor by which we can _initially_ keep a node in the cluster even though it lowers the density

# ensure the subfinal sort is always the same... we sort by size since the edge density is (roughly) constant
SUBFINAL_SORT="-k 1nr -k 5n"
BLANT_FILES=''

while echo "X$1" | grep '^X-' >/dev/null; do # first argument is an option
    case "$1" in
    -b) ONLY_BIGGEST=1; shift;;
    -1) ONLY_FIRST=1; shift;;
    -h) ONLY_BEST_ORBIT=1; shift;;
    -S) SQR_NORM=1; shift;;
    -w) WEIGHTED="$1"; W=w; shift;;
    -s) SAMPLE_METHOD="-s$2"; shift 2;;
    -s*) SAMPLE_METHOD="$1"; shift;;
    -D) BLANT_FILES="$2"; shift 2;;
    -T) LOCAL_TMP="$2"; shift 2;;
    -A) TURAN=true; shift;;
    -C) COMPLEMENT=-C; shift;;
    -v) VERBOSE=1; QUIET=''; shift;;
    -nm) PRINT_MEMBERS=0; shift;;
    -n) PRECISION="-n $2"; shift 2;;
    -n[0-9]*) PRECISION="$1"; shift;;
    -p[0-9]*) PRECISION="$1"L; shift;;
    -r) RANDOM_SEED="-r $2"; shift 2;;
    -emt) minEdgeMean="$2"; shift 2;;
    -r[0-9]*) RANDOM_SEED="$1"; shift 1;; # allow seed to be same or separate argument
    -m) minClusArg="$2"; shift 2;;
    -M) maxClus="$2"; shift 2;;
    -m[0-9]*) minClusArg="`echo $1|sed 's/^-m//'`"; shift 1;;
    -M[0-9]*) maxClus="`echo $1|sed 's/^-m//'`"; shift 1;;
    -o) OVERLAP="$2"; shift 2;;
    -o[0-9]*) OVERLAP="`echo $1|sed 's/^-o//'`"; shift 1;;
    -P[0-9]*) EDGE_PREDICT="`echo $1|sed 's/^-P//'`"; shift 1;;
    -*) die "unknown option '$1'";;
    esac
done

[ $# -ne 4 ] && die "expecting exactly 4 mandatory arguments, but you supplied $#"

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMPDIR=`mktemp -d ${LOCAL_TMP:-"/tmp"}/$BASENAME.XXXXXX`
#trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
echo "TMPDIR is $TMPDIR" >&2
[ "$BLANT_FILES" ] || BLANT_FILES="$TMPDIR"

BLANT_CMD="$1 ${PRECISION} $QUIET";
BLANT_EXE=`echo $BLANT_CMD | awk '{print $1}'`
Ks=(`echo $2 | newlines | sort -nr`); # sort the Ks highest to lowest so the below parallel runs start the higher values of k first
#[ `echo "${Ks[@]}" | wc -w` -eq 1 ] || die "no more multiple K's at the same time"
EDs=($3)
net=$4

[ -x $BLANT_EXE ] || die "'$BLANT_EXE' does not exist or is not an executable"
for k in "${Ks[@]}"; do
    [ "$k" -ge 3 -a "$k" -le 8 ] || die "One k is '$k' but must be between 3 and 8"
done

[ -f "$net" ] || die "network '$net' does not exist"
case "$net" in
*.el|*.eln) [ "$WEIGHTED" = "" ] || die "non-weighted edgelist given with -w option";;
*.elw) [ "X$WEIGHTED" = "X-w" ] || die "weighted edgelist given without -w option";;
*) die "network '$net' must be an edgeList file ending in .el";;
esac

if [ "$COMPLEMENT" != "" ]; then
    net4awk=$TMPDIR/net4awk.el$W
    #SAMPLE_METHOD=-sNBE # only method that works for now
    hawk 'NF==1{ASSERT(FNR==1, "single column only allowed on first line"); nextline}
	  NF==2{++D[$1];++D[$2];++edge[$1][$2];++edge[$2][$1]}
	  NF==3{D[$1]+=$3;D[$2]+=$3;edge[$1][$2]=$3;edge[$2][$1]=$3}
	END{
	    for(u in edge) for(v in edge) if(u<v)
		if(!(v in edge[u])) {
		    ASSERT(!(u in edge[v]),"oops");
		    if(NF==2) printf "%s\t%s\n", MIN(u,v),MAX(u,v)
		    if(NF==3) printf "%s\t%s\t%g\n", MIN(u,v),MAX(u,v),edge[u][v]
		}
	}' $net | sort -u > $net4awk
else
    net4awk=$net
fi
netCounts=`ls $net | sed 's/\.el/.count.el/'` # see if a count.el file already exists
if [ -f $netCounts -a -r $netCounts ]; then
    netCounts=`hawk '{ASSERT(FNR==1,"expecting only 1 line");n=$1;m=$2;printf "%d\t%d\t%g\n",n,m,m/choose(n,2)}' $netCounts`
else
    netCounts=`hawk 'FNR==1&&NF==1{n1=1;nextline}{++seen[$1];++seen[$2]}END{n=length(seen);m=NR-n1;printf "%d\t%d\t%g\n", n, m, m/choose(n,2)}' $net4awk`
fi
numNodes=`echo "$netCounts" | cut -f1`
numEdges=`echo "$netCounts" | cut -f2`
meanED=`echo "$netCounts" | cut -f3`

# echo "Don't forget: edgeDensity threshold for GDVs should fluctuate randomly as the cluster is built, with threshold something like ED+N(ED,ED/2), so that the expected density is still ED, but we're allowed to add nodes both above and below the threshold, at random" >&2

# Pre-run the BLANTs for each k, and re-use the files for each edge density.
BLANT_EXIT_CODE=0
if [ "$BLANT_FILES" = "$TMPDIR" ]; then
    for k in "${Ks[@]}"; do
	#minEdges=`hawk 'BEGIN{edC='$EDGE_DENSITY_THRESHOLD'*choose('$k',2);rounded_edC=int(edC); if(rounded_edC < edC){rounded_edC++;} print rounded_edC}'`
	# Use MCMC because it gives asymptotically correct concentrations *internally*, and that's what we're using now.
	# DO NOT USE -mi since it will NOT output duplicates, thus messing up the "true" graphlet frequencies/concentrations
	# Possible values: MCMC NBE EBE RES
	[ -f canon_maps/canon_list$k.txt ] || continue
	ABSOLUTE_CLIQUE_COUNT=""
	CMD=""
	$TURAN && CMD="./scripts/absolute-clique-count.sh $k $net > $TMPDIR/ACC 2>/dev/null &&"'ABSOLUTE_CLIQUE_COUNT="-A `cat $TMPDIR/ACC`";'
	CMD="$CMD $BLANT_CMD $WEIGHTED $ABSOLUTE_CLIQUE_COUNT $COMPLEMENT $RANDOM_SEED -k$k $SAMPLE_METHOD -mc$COMMUNITY_MODE $net"
	eval $CMD > $TMPDIR/blant$k.out &
    done
fi

RANDOM_SEED=`echo $RANDOM_SEED | sed 's/^-r *//'` # the "-r" is needed in CMD above but remove it for awk below.

# First argument defines the lower bound on the fractional overlap with a previous cluster beyond which a cluster
# will be REMOVED. ie., any cluster will be removed if it overlaps with a previous one by at least this fraction of its
# nodes (set to 1 to remove only EXACT duplicates)
remove-overlap-clusters(){ BOUND="$1"; shift;
    hawk ' # post-process to remove overlapping clusters
    BEGIN{ numClus=0 }
    {
	delete S;
	numNodes=$1
	edgeHits=$2;
	edgeWgts=$3;
	k=$4;
	for(i=5;i<=NF;i++) ++S[$i]
	WARN(length(S)==numNodes,"mismatch(1) in numNodes and length(S)");
	add=1;
	for(i=1;i<=numClus;i++) {
	    same=0;
	    for(u in S) if(u in cluster[i])++same;
	    if(same >= length(S)*'$BOUND'){ add=0; break;} # we check ONLY $BOUND fraction of the NEW cluster
	}
	if(add) {
	    ++numClus; edges[numClus]=edgeHits; edgeSum[numClus]=edgeWgts; kk[numClus]=k;
	    for(u in S) ++cluster[numClus][u];
	}
    }
    END{
	PROCINFO["sorted_in"]="@ind_num_asc"; # print nodes in numerical ascending order
	#PROCINFO["sorted_in"]="@ind_str_asc"; # print nodes in string ascending order
	for(i=1;i<=numClus;i++) {
	    printf "%d %d %g %d", length(cluster[i]), edges[i], edgeSum[i], kk[i]
	    for(u in cluster[i]) printf " %s", u
	    print ""
	}
    }'
}

# build actual communities; lots of RAM and CPU, eg 300GB and 24h for Skinnider 3.48M dream03
# Expects 2 arguments: edgeList.el, and "-", which is blant output sorted by count.
build-clusters(){ k=$1; shift;
    [ 3 -le "$k" -a "$k" -le 8 ] || die "build-clusters: k must be in [3,8]"
    hawk 'BEGIN{if("'$RANDOM_SEED'") srand('$RANDOM_SEED');else Srand();OFS="\t"; ID=0; edgePredict='"$EDGE_PREDICT"';
	M=1*'$minClusArg'; if(M>=1) minClus=M;
	else {
	    m='$numEdges'; n='$numNodes'; eps=m/choose(n,2);
	    if(M>0)pVal=M;else pVal=1/choose(n,2)^2; # heuristic to make pVal small enough to get significant clusters
	    for(e=1;e<m;e++)if(eps^e<pVal) break;
	    for(minClus=2;minClus<n;minClus++) if('$edgeDensity'*choose(minClus,2)>=e) break;
	    if('$VERBOSE')printf "minClus %d (%d edges = pVal %g < %g for ED %g)\n",minClus,e,eps^e,pVal,eps >"/dev/stderr"
	}
	ASSERT(minClus>=1,"minClus "minClus" must be >=1");
    }
    ARGIND==1 { # get the edge list [and weights if present]
	if(NF==1){
	    ASSERT(FNR==1, "single column only allowed on first line");
	    ASSERT(1*$1==$1,"first line has one column but must be an int");
	    Note("Got n="$1" from first line");
	    nextline; # skip the number of nodes, if present
	}
	else if(NF==2) weight=1;
	else if(NF==3) weight=$3;
	else ASSERT(0, "expecting either 2 or 3 columns but have "NF);
	degree[$1]+=weight;degree[$2]+=weight;edge[$1][$2]=edge[$2][$1]=weight
    }
    ARGIND==2 {
	if('$ONLY_FIRST' && ($2 in count)) next;
	item=$3; count[$2]+=$1; node[FNR]=$2; line[$2]=FNR;
	for(i=4; i<=NF; i++) if(edge[$2][$i] > '$minEdgeMean'/2) # heuristic: avoid too-weak edges
	    graphletNeighbors[$2][$i]=1;
	    graphletNeighbors[$i][$2]=1; # do we really need both?
    }
    function EdgesIntoS(v,       edgeHits,u) { # count edges from v into S, not including v (which can be is S)
	edgeHits=0;
	for(u in S) if(u!=v && (v in edge[u])){
	    ASSERT(edge[v][u],"edge error");
	    ++edgeHits;
	}
	return edgeHits;
    }
    function WgtEdgesIntoS(v,       edgeWgts,u) { # WEIGHT of edges from v into S, not including v (which can be is S)
	edgeWgts=0;
	for(u in S) if(u!=v && (v in edge[u])){
	    ASSERT(edge[v][u],"edge error");
	    edgeWgts+=edge[u][v];
	}
	return edgeWgts;
    }
    function highRelClusCount(u, v) { # Heuristic
	if(!(u in count) || count[u]==0) return 1;
	if(!(v in count) || count[v]==0) return 0;
	if (v in count) return count[v]/count[u]>=sqrt('$edgeDensity')*0.7; # values in [0.5,0.9] work well
	else return 1/count[u]>=0.5;
    }
    function AppendNeighbors(u,origin,    v,oldOrder, edgesIntoS) {
	oldOrder=PROCINFO["sorted_in"];
	if(u in graphletNeighbors) {
	    # The following three lines will append the neighbors in order of most edges back into S
	    # It seems reasonable, but is more expensive seems to have no significant difference. :-(
	    #for (v in graphletNeighbors[u]) edgesIntoS[v] = EdgesIntoS(v);
	    #PROCINFO["sorted_in"] = "@val_num_desc"; # for loop through edgesIntoS, largest first
	    #for (v in edgesIntoS)
	    # The default is to append them in random order... no significant difference on small networks
	    PROCINFO["sorted_in"]="randsort";
	    for (v in graphletNeighbors[u]) {
		if(!(v in visitedQ) && (!(v in line) || line[v] > line[origin]) && highRelClusCount(u, v)) {
		    QueueAdd("Q", v);
		    visitedQ[v]=1;
		}
	    }
	    PROCINFO["sorted_in"]=oldOrder;
	}
    }
    END{n=length(degree); # number of nodes in the input network
	if('$edgeDensity'==1) density_leeway=1;
	else density_leeway = '$DENSITY_LEEWAY';
	# make cluster and started empty sets; started is a list of nodes we should NOT start a new BFS on
	delete started; started[0]=1; delete started[0];
	numClus=0;
	QueueAlloc("Q");
	for(start=1; start<=int(FNR*'$edgeDensity'); start++) { # look for a cluster starting on line "start".
	    origin=node[start];
	    if(started[origin]) continue;
	    #printf "line %d origin %s;", start, origin > "/dev/stderr"
	    started[origin]=1;
	    MakeEmptySet(S);
	    MakeEmptySet(visitedQ);
	    misses=0; # how many nodes have been skipped because they did not work?
	    if(QueueLength("Q")>0){QueueDelloc("Q");QueueAlloc("Q");} #Ensure the queue is empty
	    QueueAdd("Q", origin);
	    visitedQ[origin] = 1;
	    wgtEdgeCount = edgeCount = 0;
	    while(QueueLength("Q")>0){
		u = QueueNext("Q");
		ASSERT(!(u in S),"u in S error");
		newEdgeHits = EdgesIntoS(u);
		wgtEdgeHits = WgtEdgesIntoS(u);
		edgeCount += newEdgeHits;
		wgtEdgeCount += wgtEdgeHits;
		S[u]=1;
		Slen = length(S);
		#printf "start %d u %s newEdgeHits %d edgeCount %d Slen %d\n",start,u,newEdgeHits,edgeCount,Slen>"/dev/stderr"
		ASSERT(Slen > 0, "Slen must be > 0 but is "Slen);
		if(Slen >= '$maxClus') break;
		# CAREFUL: calling InducedEdges(S) every QueueNext() is VERY expensive; uncommenting the line below
		# slows the program by more than 100x (NOT an exaggeration)
		# WARN(edgeCount == InducedEdges(edge,S),"Slen "Slen" edgeCount "edgeCount" Induced(S) "InducedEdges(edge,S));
		Sorder[Slen]=u; # no need to delete this element if u fails because Slen will go down by 1
		if (Slen==1)
		    AppendNeighbors(u, origin);
		else {
		    maxEdges = choose(Slen,2);
		    if(edgeCount/maxEdges < '$edgeDensity' * density_leeway) {
			delete S[u]; # u drops the edge density too low, so nuke it
			if(edgePredict) {
			    if((edgeCount+edgePredict)/maxEdges >= '$edgeDensity') {
				# just ONE edge would put us over the threshold; predict it should exist
				numPredict=0;
				for(uu in S)if(!(u in edge[uu]))
				    #printf "predictEdge-%d\t%s\t%s\n",++numPredict, u,uu > "/dev/stderr"
				WARN(numPredict<=edgePredict, "should only get "edgePredict" new edge predictions but got "numPredict);
			    }
			}
			edgeCount -= newEdgeHits;
			wgtEdgeCount -= wgtEdgeHits;
			if(++misses > MAX(Slen, n/100)) break; # 1% of number of nodes is a heuristic...
			# keep going until count decreases significantly; remove duplicate clusters in next awk
		    } else {
			misses=0;
			AppendNeighbors(u, orign);
		    }
		}
	    }
	    #printf " |S|=%d edgeMean %g", length(S), wgtEdgeCount/(edgeCount?edgeCount:1) > "/dev/stderr"
	    if(QueueLength("Q")==0 && length(S) > 1) {
		# post-process to remove nodes with too low in-cluster degree---quite relevant for lower density
		# communities where a node may be added because it does not
		# reduce the *mean* degree of the cluster, but it does not really have sufficiently strong
		# connections to the actual members of the community. We use two criteria:
		# 1) the in-cluster degree is more than 3 sigma below the mean, or
		# 2) the in-cluster degree is less than 1/3 the mode of the in-cluster degree distribution.
		# The latter was added in response to our performance on the LFR graphs, but it does not appear
		# to hurt performance anywhere else.
		StatReset(""); delete degFreq;
		tmpEdge = InducedEdges(edge,S, degreeInS); # this call populates degreeInS
		ASSERT(tmpEdge == edgeCount, "edgeCount "edgeCount" disagrees(1) with InducedEdges of "tmpEdge);
		for(v in S) { StatAddSample("", degreeInS[v]); ++degFreq[degreeInS[v]];}
		maxFreq=degMode=0;
		for(d=StatMax("");d>=StatMin("");d--)
		    if(d in degFreq && degFreq[d] >= maxFreq) { # use >= to extract SMALLEST mode
			maxFreq=degFreq[d]; degMode=d
		}
		#printf " deg mean %g stdDev %g maxFreq %d degMode %d; pruning...", StatMean(""), StatStdDev(""), maxFreq, degMode > "/dev/stderr"
		PROCINFO["sorted_in"] = "@val_num_asc"; # for loop through in-degrees, smallest first
		for(v in degreeInS) {
		    if(degreeInS[v] < StatMean("") - 3*StatStdDev("") || degreeInS[v] < degMode/3) {
			#printf " %s(%d)", v, degreeInS[v] > "/dev/stderr";
			delete S[v];
			edgeCount -= EdgesIntoS(v);
			wgtEdgeCount -= WgtEdgesIntoS(v);
		    }
		}
		#printf " |S|=%d", _statN[""] > "/dev/stderr"
	    }
	    tmpEdge = InducedEdges(edge,S, degreeInS); # this call populates degreeInS
	    ASSERT(tmpEdge == edgeCount, "edgeCount "edgeCount" disagrees(1) with InducedEdges of "tmpEdge);
	    PROCINFO["sorted_in"] = "@val_num_asc"; # for loop through in-degrees, smallest first
	    for(v in degreeInS) {
		if(length(S)<2) break;
		if(edgeCount / choose(length(S),2) >= '$edgeDensity') break; # break once ED is above threshsold
		delete S[v];
		edgeCount -= EdgesIntoS(v);
		wgtEdgeCount -= WgtEdgesIntoS(v);
		tmpEdge = InducedEdges(edge,S); #avoid changing the degreeInS array
		ASSERT(tmpEdge == edgeCount, "edgeCount "edgeCount" disagrees(1) with InducedEdges of "tmpEdge);
	    }
	    if(length(S)>1) {
		ed = edgeCount / choose(length(S),2);
		WARN(ed >= '$edgeDensity',"hmmm, ed is "ed);
	    }
	    #printf " |S|=%d", _statN[""] > "/dev/stderr"
	    Slen=length(S);
	    if(Slen>=minClus && wgtEdgeCount/edgeCount>='$minEdgeMean') {
		maxEdges=choose(Slen,2);
		#print " ACCEPTED" > "/dev/stderr";
		++numClus; printf "%d %d %g '$k'", Slen, edgeCount, wgtEdgeCount
		StatReset("");
		tmpEdge=InducedEdges(edge,S, degreeInS);
		ASSERT(tmpEdge == edgeCount, "edgeCount "edgeCount" disagrees(2) with InducedEdges of "tmpEdge);
		for(u in S) {printf " %s", u; StatAddSample("", degreeInS[u]);}
		#printf "final |S|=%d mean %g stdDev %g min %d max %d:\n", _statN[""], StatMean(""), StatStdDev(""), StatMin(""), StatMax("") > "/dev/stderr"
		print "";
		if('$ONLY_BIGGEST') exit;
		# now mark as "started" the first half of S that was filled... or more generally, some fraction.
		# We choose the top (1-OVERLAP) fraction because OVERLAP is meant to be less stringent as it
		# approaches 1, which in the case of THIS loop means that if OVERLAP is close to 1, we want
		# to eliminate a SMALLER proportion (ie., allow more) of future BFS starts, so if OVERLAP=1,
		# we eliminate NOTHING in this loop.
		for(i=1;i<Slen*(1-'$OVERLAP');i++) started[Sorder[i]]=1;
	    }
	    #else print " REJECTED" > "/dev/stderr";
	}
    }' "$@"
}

# check if there's a compiled executable of the appropriate name, otherwise fall back on hawk version
use-compiled() {
    EXE="$1"; shift;
    if [ -x "$EXE" ]; then # use the compiled executable
	./"$EXE" "$@"
    else # fall back to the AWK version
	"$EXE" "$@"
    fi
}

for k in "${Ks[@]}"; do
    wait; (( BLANT_EXIT_CODE += $? ))
done

for edgeDensity in "${EDs[@]}"; do
    for k in "${Ks[@]}"; do
	hawk ' # This first awk just sorts the BLANT output; it uses minimal RAM and CPU
	    BEGIN{ edC='$edgeDensity'*choose('$k',2); onlyBestOrbit='$ONLY_BEST_ORBIT';
		cMode="'$COMMUNITY_MODE'"; if(cMode=="") cMode=="g"; # graphlet uses FAR less RAM
		ASSERT(cMode=="g" || cMode=="o", "COMMUNITY_MODE must be o or g, not "cMode);
		rounded_edC=int(edC); if(rounded_edC == edC) --rounded_edC; # allow one edge less than eps
		minEdges=MAX(rounded_edC, ('$k')-1)
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
		# save the neighbors of those clusters that have high edge density for BFS
		for(j=4;j<=NF;j++) ++graphletNeighbors[$1][item][$j];
	    }
	    END {
		for(u in Kc) for(item in Kc[u]) if(item) { # do not use item==0
		    ORS=" " # Output Record Separator = space, so all the following goes on one line
		    if('$SQR_NORM') print Kc[u][item]^2/T[u], u, item
		    else            print Kc[u][item], u, item # print the cluster membership count, and the node
		    for(v in graphletNeighbors[u][item]) print v
		    ORS="\n"; print ""; # finally print the newline
		}
	    }' canon_maps/canon_list$k.txt canon_maps/orbit_map$k.txt $BLANT_FILES/blant$k.out |
	sort -T $TMPDIR -gr | # > $BLANT_FILES/clus$k.sorted  # sorted [weighted] cluster-counts of all the nodes, largest-to-smallest
	if [ -x build-clusters ]; then # use the compiled executable
	    # Help the C program by telling it how many neighbors there are since scanf() can't distinguish lines
	    # arguments: edgeDensity minClus edgeList.el blant-output
	    awk '{printf "%s %s %s %d", $1,$2,$3,NF-3; for(i=4;i<=NF;i++) printf " %s",$i; print ""}' | tee /tmp/x |
		./build-clusters $k $edgeDensity $minClusArg "$net4awk" -
	else # fall back to the AWK version
	    build-clusters "$k" "$net4awk" -
	fi |
	sort -T $TMPDIR -nr | # sort by node count, largest to smallest
	use-compiled remove-overlap-clusters 1 "$net4awk" | 
	sort -T $TMPDIR $SUBFINAL_SORT > $TMPDIR/subfinal$k$edgeDensity.out & # sort by #nodes, then by first node in list
    done
    for k in "${Ks[@]}"; do
	    wait; (( BLANT_EXIT_CODE += $? ))
    done
done

# we can use the --merge option because the subfinal files are already sorted in SUBFINAL_SORT order
sort -T $TMPDIR --merge $SUBFINAL_SORT $TMPDIR/subfinal*.out |
    use-compiled remove-overlap-clusters "$OVERLAP" "$net4awk" | 
	hawk ' # now prettyprint
	BEGIN{ PROCINFO["sorted_in"]="@ind_num_asc"; # print nodes in numerical ascending order
	    #PROCINFO["sorted_in"]="@ind_str_asc"; # print nodes in string ascending order
	}
	{ nodes=$1; edges=$2; edgeSum=$3; k=$4;
	    maxEdges=choose(nodes,2);
	    ASSERT(edgeSum/edges>='$minEdgeMean', "oops, should not get this far with an edgeMean of "edgeSum/edges);
	    printf "%d nodes, %d/%d edges, %g%% density (clusterWeight %g) from k%d", nodes, edges, maxEdges,
		100*edges/maxEdges, edgeSum, k
	    if('$PRINT_MEMBERS') {
		printf ", nodeSet {"
		for(c=5;c<=NF;c++) printf " %s", $c
		printf " }"
	    }
	    print ""
	}' |
    sort -T $TMPDIR -k 1nr -k 8gr -k 10 -k 12 # sort order: nodeCount, clusterWeight, k, nodeSet
#set -x
exit $BLANT_EXIT_CODE
