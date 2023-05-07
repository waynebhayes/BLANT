#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
measure="EDN"
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME M E t stopT measure network.el [blant.out files, or none to read stdin]
PURPOSE: run blant-clusters for E edge densities in network.el, then generate a similarity graph for it.
    M sample multiplier for blant clusters.
    E number of edge densities. It will start in 1/E and increment 1/E each step. 
    t [0,1] similarity threshold of the output communities. Communities in which the percentage of neighbors is
    t will not be retrieved.
    stopT difference in EDN increase for which it is not worth to continue expanding. 
    measure to optimize [EDN,OMOD]
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


M=$1
E=$2
t=$3
stopT=$4
measure=$5;
net=$6;

PARALLEL=${PARALLEL:-"/bin/bash"}
START_AT=${START_AT:-"0.1"}

[ -f "$net" ] || die "network '$net' does not exist"
[[ "$measure" == "EDN" || "$measure" == "OMOD" ]] || die "Measure $measure not in the list"

numNodes=`newlines < $net | sort -u | wc -l`
stepSize=$(hawk 'BEGIN{print (1-'$START_AT')/('$E'-1)}')

commands=""
for edgeDensity in $(seq -f "%.4f" $START_AT $stepSize 1.0) ; do
    commands+="./scripts/blant-clusters.sh ./blant $M '3 4 5 6 7' '$edgeDensity' $t $net > $TMPDIR/blant-c$edgeDensity.out;"
    #Printing progress so you know speed and that it is not stuck
    commands+="i=\$(find $TMPDIR -name blant-c*.out -type f -not -empty | wc -l);"
    commands+="i=\`expr 100 \"*\" \$i \`;"
    commands+="percentage=\`expr \$i / $E \`;";
    commands+="printf 'Done runninng blant-c for d=$edgeDensity[' >&2;"
    commands+="for ((j=0; j<\$percentage; j+=2)); do printf '#' >&2; done;"
    commands+="for ((j=\$percentage; j<100; j+=2)); do printf ' '>&2; done;"
    commands+="printf \"] \$percentage %%\r\">&2;\n"  
done

echo -e $commands  | $PARALLEL >&2 &
