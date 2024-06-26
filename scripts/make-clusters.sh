#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NUL=" "; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [-w] mu Ks network.el[w]
PURPOSE: call blant-clusters.sh with mu and Ks at 1% edge density increments from 100% down to 1%"

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
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

WEIGHT=''
W=''
m=''
while true; do
    case "$1" in
    -w) WEIGHT="-w"; W=w; shift;;
    -m) m="-m $2"; shift 2;;
    -*) die "unknown option '$1'";;
    *) break;;
    esac
done

[ $# -eq 3 ] || die "expecting 3 arguments"

mu="$1"; Ks="$2"; net="$3"; shift 3

echo "$mu" | grep '^[0-9][0-9]*$' >/dev/null || die "first argument must be an integer"
for k in $Ks; do echo "$k" | grep '^[3-8]$' >/dev/null || die "Ks in second argument must be from 3 to 8"; done
[ -f "$net" ] || die "no such network file '$net'"

high=1.00
low=0.99
./scripts/blant-clusters.sh $m -o 0 $WEIGHT ./blant $mu "$Ks" $high "$net" | tee $TMPDIR/$high.out |
    sed 's/.*://' | gawk 'ARGIND==1{for(i=1;i<=NF;i++)++taken[$i]}ARGIND==2 && !($1 in taken) && !($2 in taken){print}' - "$net" > $TMPDIR/net-$high.el$W
cat $TMPDIR/$high.out
count.el $TMPDIR/net-$high.el$W

for d in `seq -f "%02g" 99 -7 1`; do
    low=0.$d
    ./scripts/blant-clusters.sh $m -o 0 $WEIGHT ./blant $mu "$Ks" $low $TMPDIR/net-$high.el$W | tee $TMPDIR/$low.out |
	sed 's/.*://' | gawk 'ARGIND==1{for(i=1;i<=NF;i++)++taken[$i]}ARGIND==2 && !($1 in taken) && !($2 in taken){print}' - $TMPDIR/net-$high.el$W > $TMPDIR/net-$low.el$W
    cat $TMPDIR/$low.out
    count.el $TMPDIR/net-$low.el$W
    high=$low
done
