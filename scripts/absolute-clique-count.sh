#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NUL=" "; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME k graph.el
PURPOSE: use Sweta Jain's clique-counting code to get an estimate of the number of k-cliques across entire graph.el
    We assume that her code is in the directory src/cliquecounting"

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

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

# Create an "edges" file for Sweta Jain's code, which is just an edgelist preceded by a line with node and edge count.
[ $# -eq 2 ] || die "expecting exactly 2 arguments"
[ $1 -ge 3 -a $1 -le 8 ] || die "first argument must be an integer from 3 to 8 but is '$1'"
[ -f "$2" ] || die "cannot find file '$2'"
[ -d src/cliquecounting/tests ] || die "can't find Sweta Jain's code directory"
[ -x src/cliquecounting/tests/test_cliques_turan_shadow ] || die "can't find Jain's executable"
k=$1; shift;

# Unfortunately her code expects files to be in a particular place:
mkdir -p src/cliquecounting/graphs src/cliquecounting/results/cliques || die "failed to make Turan directories"
F=src/cliquecounting/Escape/nCr.txt; [ -f $F ] || unxz <$F.xz >$F || die "failed to create nCr.txt"
(
    B="`basename "$1"`"
    # Now ensure node names are integers
    awk '!id[$1]{name[++n]=$1;id[$1]=n}!id[$2]{name[++n]=$2;id[$2]=n}{m++;edge[$1][$2]=1}
	END{print n,m; for(u in edge)for(v in edge[u])print id[u]-1,id[v]-1}' "$1" > src/cliquecounting/graphs/"$B"
    r=1; n=500000; cd src/cliquecounting/tests; ./test_cliques_turan_shadow -m $k -M $k -r $r -n $n -i "$B" 2>&1 | grep -i 'number of.*cliques' | awk 'BEGIN{max=0}$NF>max{max=$NF}END{printf "%d\n", max}'
)
