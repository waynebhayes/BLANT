#!/bin/bash
################## SKELETON: DO NOT TOUCH THESE 2 LINES
BASENAME=`basename "$0" .sh`; TAB='	'; NL='
'
measure="EDN"
#################### ADD YOUR USAGE MESSAGE HERE, and the rest of your code after END OF SKELETON ##################
USAGE="USAGE: $BASENAME t stopT measure network.el [blant.out files, or none to read stdin]
PURPOSE: run blant-clusters for E edge densities in network.el, then generate a similarity graph for it.
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


t=$1
stopT=$2
measure=$3;
net=$4;
shift 4

[ -f "$net" ] || die "network '$net' does not exist"
[[ "$measure" == "EDN" || "$measure" == "OMOD" ]] || die "Measure $measure not in the list"

numNodes=`newlines < $net | sort -u | wc -l`

sort -k 1nr -k 3nr -k 11n "$@" | # sed 's/[^0-9.]/ /g' | # for C version, remove everything not a number
hawk 'BEGIN{ numCliques=0 } # post-process to remove duplicates
			{
			delete S; 
			numNodes=$1;
			edgeHits=$3;
			maxEdges=$5;
			k=$9;
			for(i=11;i<=NF;i++) ++S[$i]
			ASSERT(length(S)==numNodes,"mismatch in numNodes and length(S)");
			add=1;
			for(i=1;i<=numCliques;i++) {
					same=0;
					for(u in S) if(u in cluster[i]) ++same;
					if(same > length(cluster[i])*'$t'){ add=0; break;}
			}
			if(numCliques==0 || add==1) {
					++numCliques; 
					printf "%d nodes, %d of %d edges from k %d (%g%%):",
				length(S), edgeHits, maxEdges, k, 100*edgeHits/maxEdges
					for(u in S) {cluster[numCliques][u]=1; printf " %s", u}
					print ""
			}
			}' | tee $TMPDIR/blant.out |
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
    }' | 
hawk 'BEGIN{ Q=0; delete s;Srand(); measure="'$measure'"}
      ARGIND==1{++degree[$1];++degree[$2];A[$1][$2]=A[$2][$1]=1}
      ARGIND==2{neighbors[$1][$2]=neighbors[$2][$1]=$3}                                    #s number of times node u appears
      ARGIND==3{
        for(i=11;i<=NF;i++){comm[FNR][$i]=1}; 
        edges[FNR]=$3; 
        edgeDensity[FNR]=$3/$5;
        nc[FNR]=length(comm[FNR])
        for (u in comm[FNR]){
          kinu=0;
          for (v in A[u]){if(v in comm[FNR]) kinu++}
          kin[FNR][u]=kinu
        }
      }
      function scoreOfNodeInCommunity(c,u,ms){
        if(measure=="OMOD"){
          return ( (kin[c][u] - ( degree[u]-kin[c][u] ) ) / degree[u] ) * ( 1 / ms ) * edgeDensity[c] * ( 1 / nc[c] )
        } else{
          return ( 1 / ms ) * edgeDensity[c]
        }
      }
      function communityScore(c,       CS){
        if (c in score) return score[c];
        CS=0;
        for (u in comm[c]){
          CS+=scoreOfNodeInCommunity(c, u, s[u]+1)
        }
        score[c]=CS
        return CS
      }
      function potentialScore(c,       P){
        P=Q
        for(u in comm[c]){
            if(!(u in s)) continue
            for (c2 in finalComm){
              if(u in comm[c2]){
                P=P-scoreOfNodeInCommunity(c2, u, s[u])+scoreOfNodeInCommunity(c2, u, s[u] + 1)   
              }
            }
        }
        P+=communityScore(c)
        return P
      }
      function addToResult(c){
        P=potentialScore(c)
        diff=P-Q
        if (diff < '$stopT') return;
        Q=P
        for (u in comm[c]){s[u]++}
        finalComm[c]=1
      }
      function markDFS(c){
        visitedDFS[c]=1;connectedSet[c]=1;
        if (!(c in neighbors)) return;
        for(n in neighbors[c]){
          if (n in visitedDFS) continue;
          markDFS(n);
        }
      }
      function getRandom(set){
        i=int(rand()*length(set))
        count=0
        for (c in set){
          if (count==i) return c
          count++
        }
      }
      function getBiggest(set){
        size=0; b=-1
        for (c in set){
          if (nc[c]>size){size=nc[c];b=c}
        }
        return b
      }
      function expand(c){
        if (!(c in neighbors)) return;
        for (n in neighbors[c]){
          if (n in visitedComm) continue;
          P=potentialScore(n);
          diff=P-Q
          if ( (diff > '$stopT') ){
              PQpush("PQ",P,n)
          }
          visitedComm[n]=1  
        }
      }
      END{
        N=length(degree); K=length(comm)
        delete visitedComm; delete visitedDFS;
        for (c in comm){
          if(c in visitedDFS) continue
          delete connectedSet;
          markDFS(c)
          cbest=getRandom(connectedSet)
          PQpush("PQ", communityScore(cbest),cbest)
          visitedComm[cbest]=1;
        }
        while (PQlength("PQ")>0){
          delete score;
          c=PQpop("PQ")
          addToResult(c)
          expand(c)          
        }
        if(measure=="OMOD"){
          printf "Qov=%s",Q/length(finalComm)
        } else{
          printf "EDN=%s",Q
        }
        print ""
        for (c in finalComm){
            maxEdges=choose(nc[c],2)
            printf "%d nodes, %d of %d edges (%g%%):", nc[c], edges[c], maxEdges, 100*edgeDensity[c]
            for(u in comm[c]) {printf " %s", u}
            print ""
        }
      }' $net - $TMPDIR/blant.out | sort -k 1nr -k 11n
