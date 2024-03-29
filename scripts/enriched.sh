#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NUL=" "; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME [-k K] [-p] [-G network.el ] searchRegexp eraseRegexp pTresh taxID GENE2GO clusters-file(s)
PURPOSE: given a gene2go file (note we do NOT perform hiararchy expansion yet), and a file containing clusters of proteins
    (one clutser per line) output clusters that are enriched with p-value less than pThresh according to the Hypergeometric
    test, and optionally [using the -p option] print any proteins that do NOT have the GO term on a 'prediction' line.
ARGUMENTS
    searchRegexp - perform enrichment analysis only on lines matching this regexp, eg 'nodes,'
    eraseRegexp - the regexp matching the part of the line that should be erased, eg '^.*: '
    pThresh - only print GO terms in a cluster if its enrichment p-value is less than this (can be an AWK expression)
    taxID - eg 9606 for human
OPTIONS
    -p - print predictions ripe for 'fgrep -f' on a newer gene2go file
    -k k-val - only include clusters coming from graphlets of size k or larger (default 3, ie. unrestricted)
    -G network.el - used optionally with -p to print a GO annotation prediction only if some fraction of the predicted
	node's in-cluster neighbors are annotated [use with caution, as it doesn't appear to improve anything and may
	be counter-productive since the whole point is we're using the COMMUNITY to make the prediction, not individual
	neighbors]"

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

PREDICT=0
USE_GRAPH=0
k=3
G=''
while true; do
    case "$1" in
    -p) PREDICT=1; shift;;
    -G) USE_GRAPH=1; G="$2"; shift 2;;
    -k) k="$2"; shift 2;;
    -*) die "unknown option '$1'";;
    *) break;;
    esac
done

if [ "$G" != "" ]; then
    [ -f "$G" ] || die "cannot open network file '$G'"
    [ "$PREDICT" = 1 ] || die "network file '$G' is useless unless prediction is enabled with the -p option"
fi

[ $# -ge 6 ] || die "expecting at least 6 arguments after options"

sRE="$1"
eRE="$2"
pThresh="$3"
taxID="$4"
shift 4

hawk 'BEGIN{if('$USE_GRAPH') MakeEmptySet(edge)}
    '$USE_GRAPH'==1&&ARGIND==1{edge[$1][$2]=edge[$2][$1]=1}
    ARGIND=='$USE_GRAPH'+1&&$1=='$taxID'{++pGO[$2][$3];++pGOev[$2][$3][$4];++GOp[$3][$2]}
    ARGIND>'$USE_GRAPH'+1 && /'"$sRE"'/{
	print; sub("'"$eRE"'",""); delete annotGP; delete annotPG;
	for(c=1;c<=NF;c++) if($c in pGO) {p=$c;for(g in pGO[p]){++annotGP[g][p];++annotPG[p][g]}}
	for(g in annotGP){
	    s=length(annotGP[g]); n=NF; S=length(GOp[g]); N=length(pGO); pv=HyperGeomTail(s,n,S,N);
	    if(pv<('"$pThresh"')) {
		printf "\t%s (%d %d %d %d) = %g\n", g,s,n,S,N,pv
		if('$PREDICT' && s<n) { # we can predict some GO terms in this community
		    for(c=1;c<=NF;c++) if(!($c in GOp[g]) && ($c in edge)) {
			numNeigh = numAnnotNeigh = 0; # Count the number of g-annotated in-cluster neighbors
			for(d=1;d<=NF;d++) if($d!=$c && ($d in edge[$c])) {
			    ++numNeigh; if($d in GOp[g]) ++numAnnotNeigh;
			}
			ASSERT(numNeigh>0, "um... "$c" has no neighbors in the cluster???");
			if(numAnnotNeigh/numNeigh > 0.5*length(annotGP[g])/NF) # fraction of g-annot neigh > frac annot in clus
			    printf "    %d %d %d %d predict\t'$taxID'\t%s\t%s\t\n",s,n,S,N,$c,g;
		    }
		}
	    }
	}
    }' $G "$@" | awk '!seen[$0]{print}{++seen[$0]}'
