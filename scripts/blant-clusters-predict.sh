#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
EXEDIR=`dirname "$0"`; BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME bla bla bla
PURPOSE: describe purpose in words"

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

[ "$MYTMP" ] || export MYTMP="/tmp"
export TMPDIR=${TMPDIR:-`mktemp -d $MYTMP/$BASENAME.XXXXXX`}
#trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

#################### END OF SKELETON, ADD YOUR CODE BELOW THIS LINE

[ $# -eq 3 ] || die "usage: taxID density date"

tx=$1
density=$2
dateDash=`echo $3 | sed 's,/,-,'`
dateDir=`echo $3 | sed 's,-,/,'`

fgrep $dateDash $HOME/extra1/preserve/BioGRID/release-dates.tsv || die "no release for $3"
R=`fgrep $dateDash $HOME/extra1/preserve/BioGRID/release-dates.tsv | sed 's/-..$//' | tail -1 | awk '{print $1}'`

net=../SANA/networks/HSapiens-$R.el

shift 3

[ -f $net ] || die "no such network $net"

GO1=$HOME/extra1/preserve/GO/processed/Entrez/$dateDir/gene2go.NOSEQ
GO2=$HOME/extra1/preserve/GO/processed/Entrez/2023/04/gene2go.NOSEQ

[ -f "$GO1" ] || die "no $GO1"
[ -f "$GO2" ] || die "no $GO2"

./scripts/blant-clusters.sh ./blant.k4-good '3 4 5' $density $net |
    while read K line; do
    [ $K -lt 10 ] && break; echo "##################### $K $line" | sed 's/{.*}//'
    echo "$line" | sed 's/.*{//' -e 's/}//' | newlines | tee $TMPDIR/clus.txt |
	awk '{printf "'$tx'\t%d\t\n",$1}' | fgrep -f - $GO1 | cut -f2-3 | sort -u > $TMPDIR/prot-GO.txt
	cut -f2 $TMPDIR/prot-GO.txt | sort -u > $TMPDIR/GO.txt # uniq list of GO terms already in the cluster
	cut -f2 $TMPDIR/prot-GO.txt | sort | uniq -c | sort -nr | # GO terms sorted by frequency in the cluster
	awk "\$1>$K/2" |
	while read k g; do
	    awk '{printf "'$tx'\t%s\t'$g'\t\n", $1}' $TMPDIR/clus.txt | fgrep -f - $GO2 |
		fgrep -vf <(awk '{printf "\t%s\t\n", $0}' $TMPDIR/prot-GO.txt) | cut -f2-3 | sort -u | # cut out prot-GO
		tee $TMPDIR/p-g-validated.txt | cut -f1 | sort -u | tee $TMPDIR/p-validated.txt |
		hawk 'ARGIND==1&&/^'$tx'	/{++p[$2];++GOp[$3][$2]}ARGIND==2{nVal=FNR}
		    END{nPred='$K-$k'; lambda=length(GOp["'$g'"]);
			pVal=HyperGeomTail(nVal,nPred,lambda,length(p));
			if(nVal && pVal <= 1/lambda)
			    printf "'$K'-node community, '$g' annotates '$k' : %d predictions, %d validated [%d%%] (lambda %d, p=%g)\n",
				nPred, nVal,100*nVal/(nPred), lambda, pVal
			}' $GO2 -
	done
    done

