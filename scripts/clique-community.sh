#!/bin/sh
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME N k
PURPOSE: output a graph containing N cliques of size k in which every node is additionally connected to exactly one node
in every other clique. If N > k, this results in every node having more edges heading outside its clique than inside, which
violates some definitions of 'community', even though most reasonable people would agree that these cliques constitute
communties inside the larger network."

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
[ "$MYTMP" ] || MYTMP=`for i in /scratch/preserve/wayne /var/tmp/wayne /tmp/wayne; do mkdir -p $i && break; done; echo $i`
TMPDIR=`mktemp -d $MYTMP/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

[ $# -eq 2 ] || die "expecting exactly 2 arguments"

N="$1"
k="$2"

cat /dev/null | # ensure only the BEGIN statement gets executed
    awk 'BEGIN{
	N='$N';k='$k';
	for(c=0;c<N;c++) {
	    for(i=0;i<k;i++) for(j=i+1;j<k;j++){
		u=c*k+i; v=c*k+j;print u,v
	    }
	}
	for(c=0;c<N;c++) for(d=c+1;d<N;d++) for(i=0;i<k;i++){
	    u=k*c+i;v=k*d+(i+c)%k;
	    print u,v
	}
    }'
