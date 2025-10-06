#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME N k orbit-pair train.el test.el
PURPOSE: All arguments are required:
    N = how many predictions to make
    k = size of graphlets to use
    orbit-pair = just like it sounds (eg 11:11 for the L3-path of Kovacs et al)
    train.el, test.el: the training and test networks, respectively"

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

[ $# -eq 5 ] || die "wrong number of arguments"

# compile the file src/blant-predict-release.c using libwayne libraries:
TOP="$1" # number of predictions to take
k=$2
OP="$3"
TRAIN="$4"
TEST="$5"

LIBWAYNE_HOME=`cd libwayne && /bin/pwd`
PATH="$PATH:$LIBWAYNE_HOME/bin:$LIBWAYNE_HOME/../scripts"
export LIBWAYNE_HOME PATH

echo "Note: BLANT will print status messages to stderr about its use of batch sampling to estimate orbit-pair counts to about 1 digit of precision; then it will print the corroborated predictions out of the top $TOP predictions as you requested; pipe thetoutput to 'wc -l' to count them" | fmt >&2 

# Now run predictions
./blant -p1L -sMCMC -k$k -mp$OP "$TRAIN" |
    sort -gr | # sort the output node-pairs by OP count, highest-to-lowest
    gawk 'BEGIN{N='$TOP'} # store how many predictions to extract
	ARGIND==1{train[$1][$2]=train[$2][$1]=1} # training edges from $TMP.train, given as first filename argument
	ARGIND==2{ test[$1][$2]= test[$2][$1]=1} # test edges from $TMP.test, given as second filename argument
	ARGIND==3{ # the third filename argument is the 3 columns from the above pipeline: score, node1, node2
	    # Note we ignore the score in $1; we look only at the 2nd and 3rd columns ($2,$3)
	    if(($2 in train) && ($3 in train[$2])) next; # skip training edges
	    ++numPred; # this line was not in the training set, and so counts as a prediction
	    if(numPred > N) exit; # exit once we have made N predictions
	    if(($2 in  test) && ($3 in  test[$2])) print; # print matching test edges
	}' "$TRAIN" "$TEST" - # the three filename arguments (dash refers to stdin, which is from the pipeline above)
