#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME network.el M E t
PURPOSE: run blant-clusters for E edge densities in network.el, then generate a similarity graph for it.
    M sample multiplier for blant clusters.
    E number of edge densities. It will start in 1/E and increment 1/E each step. 
    t [0,1] similarity threshold of the output communities. Communities in which the percentage of neighbors is
    t will not be retrieved.
    stopT difference in EDN increase for which it is not worth to continue expanding. 
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

net=$1 
M=$2
E=$3
t=$4
stopT=$5
edgeDensityStep=`hawk "BEGIN{print 1/$E}"`
edgeDensity=$edgeDensityStep

[ -f "$net" ] || die "network '$net' does not exist"

while [ $E -gt 0 ]; do
    ./scripts/blant-clusters.sh ./blant $M $net $t $edgeDensity > $TMPDIR/blant$edgeDensity.out
    edgeDensity=`hawk "BEGIN{print $edgeDensity+$edgeDensityStep}"`
    E=$((E-1))
done


sort -k 1nr -k 11n $TMPDIR/blant*.out > $TMPDIR/blant.out
numNodes=`newlines < $net | sort -u | wc -l`
hawk 'BEGIN{delete intersection}
    {
    for(i=11;i<=NF;i++){comm[FNR][$i]=1}
    for (i=1; i<FNR; i++){
        SetIntersect(intersection,comm[i], comm[FNR])
        similarity=length(intersection)
        if (similarity > 0){
            printf "%d %d %d", i, FNR, similarity
            print ""
        }
    }
    }' $TMPDIR/blant.out | 
hawk 'BEGIN{delete comm; EDN=0; delete d; delete s; delete finalComm;}     #d current density added by node u. 
      ARGIND==1{neighbors[$1][$2]=neighbors[$2][$1]=$3}  #s number of times node u appears
      ARGIND==2{for(i=11;i<=NF;i++){comm[FNR][$i]=1}; edges[FNR]=$3; edgeDensity[FNR]=$3/$5; score[FNR]=edgeDensity[FNR]*$1}
      function communityScore(c,    CS){
        if (c in score) return score[c]
        CS=0
        for (u in comm[c]){CS+=1/(s[u]+1)}
        CS*=edgeDensity[c]
        score[c]=CS
        return CS
      }
      function potentialScore(c,       P){
        P=EDN
        for(u in comm[c]){
            if (u in d) P= P - d[u] * 1/s[u] + d[u] * 1 / (s[u] + 1)
        }
        P+=communityScore(c)
        return P
      }
      function addToResult(c){
        P=potentialScore(c)
        visitedComm[c]=1
        if ((P-EDN) < '$stopT') return
        EDN=P
        for (u in comm[c]){s[u]++; d[u]+=edgeDensity[c]}
        finalComm[c]=1
      }
      function getMaximum(c,       score, best, m, n){
        score=communityScore(c)
        visitedComm[c]=1; best=c;
        if (!(c in neighbors)) return best;
        for(n in neighbors[c]){
            if (n in visitedComm) continue
            m=getMaximum(n)
            if(communityScore(m)>score) best=m
        }
        return best
      }
      function expand(c){
        if (length(s)>='$numNodes') return
        if (!(c in neighbors)) return;
        for (n in neighbors[c]){
            if (n in visitedComm) continue
            P=potentialScore(c)
            if (((P-EDN) > '$stopT') && neighbors[c][n]/MAX(length(comm[c]),length(comm[n]))  < '$t'){
                QueueAdd("Q", n)
            }
            visitedComm[n]=1
        }
      }
      END{ 
        QueueAlloc("Q");
        delete visitedComm
        for(c in comm){
            if(visitedComm[c]) continue
            b=getMaximum(c) #DFS to get maximum of every disconnected relationship subgraph
            QueueAdd("Q",b)
        }
        delete visitedComm;
        while (QueueLength("Q")>0){
            delete score;
            c=QueueNext("Q")
            addToResult(c)
            expand(c)
        }
        printf "EDN=%s",EDN
        print ""
        for (c in finalComm){
            maxEdges=choose(length(comm[c]),2)
            printf "%d nodes, %d of %d edges (%g%%):", length(comm[c]), edges[c], maxEdges, 100*edgeDensity[c]
            for(u in comm[c]) {printf " %s", u}
            print ""
        }
      }
      
      ' - $TMPDIR/blant.out
