#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [blant-clusters-output-file(s)]
PURPOSE: given the output (possibly piped) of one or more blant-clusters.sh output files, merge nearby clusters that
    overlap (even if the result is below the previously requested edge density), and spit out fewer, larger clusters
    that do NOT significantly overlap anymore."

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

[ "$MYTMP" ] || export MYTMP="/tmp"
export TMPDIR=${TMPDIR:-`mktemp -d $MYTMP/$BASENAME.XXXXXX`}
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

NODES_LOWER_BOUND=20 # ignore any and all clusters smaller than this

fgrep -h nodes, "$@" | sort -nr | # sort largest-to-smallest
    sed 's/.*{ //' -e 's/ }//' | # remove everything except the list nodes, one cluster per line
    hawk 'NF>='$NODES_LOWER_BOUND'{
	for(c=1;c<=NF;c++) ++clus[FNR][$c]; # record this cluster
	merged=0; # find out if we should merge it with a previous cluster
	    for(i=1;i<FNR;i++) {
		o=SetIntersect(res,clus[i],clus[FNR]);
		printf "cluster %d isarray %d; overlap %d\n", i, isarray(clus[i]), o;
		if(isarray(clus[i]) && (n=SetIntersect(res,clus[i],clus[FNR])) > 0) {
		    merged=i; # merge it with cluster i
		    printf "Line %d will merge with line %d since it overlaps in %d nodes\n", FNR,i,n
		    for(j=1;j<=NF;j++) ++clus[i][$j]; # merge all current nodes into cluster i
		    delete clus[FNR]; # delete this cluster since it was merged
		    next; # move on to next line
		}
	    }
	}
    END {
	for(i in clus) {
	    printf "%d nodes initially from line %d, nodeSet {", length(clus[i]), i
	    for(j in clus[i]) printf " %s", j
	    print " }"
	}
    }'
