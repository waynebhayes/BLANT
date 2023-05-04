#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
E=10
USAGE="USAGE: $BASENAME [OPTIONS] blant.exe M network.el t [E: Edge densities we are running blant-c for, default $E]
PURPOSE: use random samples of k-graphlets from BLANT in attempt to find large clusters in network.el.
    blant.exe is the name of the executable BLANT to use (usually just './blant')
    M is the mean number of times each *node* should be touched by a graphlet sample,
    so BLANT will thus take M*(n/k) total samples of k-node graphlets. 
    t similarity threshold of the output communities
    E is optional and defaults to $E.
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

WEIGHTED=0
ONLY_ONE=0
exclusive=0
SAMPLE_METHOD=-sMCMC #-sNBE
ONLY_BEST_ORBIT=0
PARALLEL=${PARALLEL:-"/bin/bash"}

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

[ $# -lt 3 ] && die "not enough arguments"
#[ $# -gt 5 ] && die "too many arguments"

BLANT=$1;
Ks=(7 6 5 4 3) #(`echo $2 | newlines | sort -nr`); # sort the Ks highest to lowest so the below parallel runs start the higher values of k first
#[ `echo "${Ks[@]}" | wc -w` -eq 1 ] || die "no more multiple K's at the same time"

sampleMultiplier=$2
net=$3;
t=$4
if [ $# -eq 5 ]; then
    E=$5;
fi

numNodes=`newlines < $net | sort -u | wc -l`

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
    CMD="$BLANT -k$k -n$n $SAMPLE_METHOD -mc $net" #-e$minEdges
    $CMD > $TMPDIR/blant$k.out &
done

for k in "${Ks[@]}"; do
    wait; (( BLANT_EXIT_CODE += $? ))
done

edgeCount=`hawk '{delete line; line[$1]=1; line[$2]=1; for (edge in line) printf "%d ",edge; print ""}' $net | sort -k 1n -k 2n | uniq | wc -l`
graphEd=`hawk 'BEGIN{print '$edgeCount'/(choose('$numNodes',2))}'`
stepSize=$(hawk 'BEGIN{print (1-'$graphEd')/('$E'-1)}')
commands=""
for edgeDensity in $(seq -f "%.4f" $graphEd $stepSize 1.0) ; do
    commands+="./scripts/community-discovery.sh $edgeDensity $ONLY_BEST_ORBIT $WEIGHTED $ONLY_ONE $TMPDIR $net $t \n"
done

echo -e $commands | $PARALLEL

sort -k 1nr -k 3nr -k 11n $TMPDIR/final*.out |
hawk 'BEGIN{ numCliques=0 } # post-process to remove duplicates
			{
			delete S; 
			numNodes=$1;
			edgeHits=$3;
			maxEdges=$5;
			k=$9;
			for(i=11;i<=NF;i++) ++S[$i]
			ASSERT(length(S)==numNodes,"mismatch in numNodes and length(S)");
			add=1;
			for(i=1;i<=numCliques;i++) {
					same=0;
					for(u in S) if(u in cluster[i]) ++same;
					if(same > length(cluster[i])*'$t'){ add=0; break;}
			}
			if(numCliques==0 || add==1) {
					++numCliques; 
					printf "%d nodes, %d of %d edges from k %d (%g%%):",
				length(S), edgeHits, maxEdges, k, 100*edgeHits/maxEdges
					for(u in S) {cluster[numCliques][u]=1; printf " %s", u}
					print ""
			}
			}'
#set -x
exit $BLANT_EXIT_CODE
