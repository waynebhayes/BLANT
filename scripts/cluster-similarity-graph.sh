#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME network.el M E t
PURPOSE: run blant-clusters for E edge densities in network.el, then generate a similarity graph for it.
    M sample multiplier for blant clusters
    E number of edge densities. It will start in 1/E and increment 1/E each step. 
    t [0,1] similarity threshold of the output communities. Communities in which the percentage of neighbors is
    t will not be retrieved
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
[ $# -lt 4 ] && die "not enough arguments"

net=$1 
M=$2
E=$3
t=$4
edgeDensityStep=`hawk "BEGIN{print 1/$E}"`
edgeDensity=$edgeDensityStep

[ -f "$net" ] || die "network '$net' does not exist"

while [ $E -gt 0 ]; do
    ./scripts/blant-clusters.sh ./blant $M $net $t $edgeDensity > $TMPDIR/blant$edgeDensity.out
    edgeDensity=`hawk "BEGIN{print $edgeDensity+$edgeDensityStep}"`
    E=$((E-1))
done


sort -k 1nr -k 11n $TMPDIR/blant*.out > $TMPDIR/blant.out
cat $TMPDIR/blant.out

hawk 'BEGIN{delete intersection}
    {
    for(i=11;i<=NF;i++){clique[FNR][$i]=1} #read clique
    for (i=1; i<FNR; i++){
        SetIntersect(intersection,clique[i], clique[FNR])
        similarity=length(intersection)
        if (similarity > 0){
            printf "%d %d %d", i, FNR, similarity
            print ""
        }
    }
    }' $TMPDIR/blant.out >&2
