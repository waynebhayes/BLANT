#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME N k [ED]
PURPOSE: output a graph containing N cliques of size k in which every node is additionally connected to exactly one node
in every other clique. If N > k, this results in every node having more edges heading outside its clique than inside, which
violates some definitions of 'community', even though most reasonable people would agree that these cliques constitute
communties inside the larger network. One can optionally provide an edge density < 1 for the 'cliques', though the script
will fail if the provided ED is too low to maintain a reasonable density of the 'cliques' inside the larger network."

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

[ $# -eq 2 -o $# -eq 3 ] || die "expecting exactly 2 arguments"

N="$1"
k="$2"
ED=1
if [ $# -eq 3 ]; then
    ED="$3"
fi

IN=`parse.awk "$ED*$N*choose($k,2)"`;
OUT=`parse.awk "choose($N,2)*$k"`;
TOT=`parse.awk "choose($k*$N,2)"`;
meanEDpc10=`parse.awk "10*int(100*($IN+$OUT)/$TOT)"`
EDpc=`parse.awk "int(100*$ED)"`
if [ $EDpc -gt $meanEDpc10 ]; then
    cat /dev/null | # ensure only the BEGIN statement gets executed
	awk 'BEGIN{
	    N='$N';k='$k';
	    for(c=0;c<N;c++) { # create the communities
		for(i=0;i<k;i++) for(j=i+1;j<k;j++) if(rand() <= '$ED') {
		    u=c*k+i; v=c*k+j; print u,v
		}
	    }
	    for(c=0;c<N;c++) for(d=c+1;d<N;d++) for(i=0;i<k;i++){ # create external edges
		u=k*c+i;v=k*d+(i+c)%k;
		print u,v
	    }
	}' > $TMPDIR/cc.el

    count.el $TMPDIR/cc.el |
	hawk '{n=$1;m=$2; meanED=m/choose(n,2);
	    if('$ED' < 2*meanED) {
		printf "ERROR: you have requested communities with edge density '$ED', which is not sufficiently higher than the mean edge density %g of the network\n", meanED > "/dev/stderr"
		exit 1;
	    }
	}' && cat $TMPDIR/cc.el
else
    parse "($IN+$OUT)/$TOT" >&2
    die "in-clique ED ($ED) is not sufficiently higher than mean ED"
fi
