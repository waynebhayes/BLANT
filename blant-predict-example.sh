#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
NET=networks/HI-union.el
USAGE="USAGE: $BASENAME N [network.el]
PURPOSE: run $BASENAME on network.el [default $NET]: split into train + test (90% & 10% of edges, respectively), and
    then count the number of correct predictions (out of N) in the test set after training."

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
which(){ echo "$PATH" | tr : "$NL" | awk '!seen[$0]{print}{++seen[$0]}' | while read d; do eval /bin/ls $d/$N; done 2>/dev/null | newlines; }
HardPath(){ if [ -h "$1" ]; then link=`cd $(dirname "$1") && /bin/ls -l $(basename "$1") | awk '/ -> /{print $NF}'`; (cd $(dirname "$1") && HardPath "$link"); else echo $(cd $(dirname "$1")&&/bin/pwd)/$(basename "$1"); fi;}

export TMPDIR="`mktemp -d /tmp/$BASENAME.XXXXXX`"
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

# compile the file src/blant-predict-release.c using libwayne libraries:
if [ $# -lt 1 -o $# -gt 2 ]; then die "must provide number of predictions as first argument"; fi
TOP="$1" # number of predictions to take

LIBWAYNE_HOME=`cd libwayne && /bin/pwd`
PATH="$PATH:$LIBWAYNE_HOME/bin:$LIBWAYNE_HOME/../scripts"
export LIBWAYNE_HOME PATH
# compile the executable
[ -f libwayne/libwayne.a ] || (cd libwayne && make all)
wgcc -O3 -Isrc -o blant-predict-release src/blant-predict-release.c

[ $# -eq 2 ] && NET="$2"
BASE_NET="`basename "$NET" .el`"

# create a random-order list of network's edges
TMP=$TMPDIR/$BASE_NET
randomizeLines < "$NET" > $TMP.el

# separate into 90% train, 10% test
FOLD=`wc -l < $TMP.el` # number of edges in network
FOLD=`expr $FOLD / 10` # take 10% for testing
head -$FOLD $TMP.el > $TMP.test
FOLD1=`expr $FOLD + 1` # training size
tail -n +$FOLD1 $TMP.el > $TMP.train

# Now run predictions
./blant-predict-release -k4 -mp6,0,0 $TMP.train | # estimate L3-paths on 4-node graphlets as in Kovacs et al.
    sort -gr | # sort the output node-pairs by L3 count, highest-to-lowest
    gawk 'BEGIN{N='$TOP'} # store how many predictions to extract
	ARGIND==1{train[$1][$2]=train[$2][$1]=1} # training edges from $TMP.train, given as first filename argument
	ARGIND==2{ test[$1][$2]= test[$2][$1]=1} # test edges from $TMP.test, given as second filename argument
	ARGIND==3{ # the third filename argument is the 3 columns from the above pipeline: score, node1, node2
	    # Note we ignore the score in $1; we look only at the 2nd and 3rd columns ($2,$3)
	    if(($2 in train) && ($3 in train[$2])) next; # skip training edges
	    ++numPred; # this line was not in the training set, and so counts as a prediction
	    if(numPred > N) exit; # exit once we have made N predictions
	    if(($2 in  test) && ($3 in  test[$2])) print; # print matching test edges
	}' $TMP.train $TMP.test - | # the three filename arguments (dash refers to stdin, which is from the pipeline above)
    wc -l # count the number of corroborated predictions
