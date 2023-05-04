#!/bin/bash

Ks=(7 6 4 5 3)
edgeDensity=$1;
ONLY_BEST_ORBIT=$2;
WEIGHTED=$3;
ONLY_ONE=$4;
TMPDIR=$5
net=$6
t=$7

for k in "${Ks[@]}"; do
    hawk 'BEGIN{ edC='$edgeDensity'*choose('$k',2); onlyBestOrbit='$ONLY_BEST_ORBIT';
        rounded_edC=int(edC); if(rounded_edC < edC) rounded_edC++;
        minEdges=rounded_edC
        }
        ARGIND==1 && FNR>1 && $2 {canonEdges[FNR-2]=$3}
        ARGIND==2 && FNR>1 && ((FNR-2) in canonEdges) {for(i=1;i<=NF;i++)orbit2canon[$i]=FNR-2; canon2orbit[FNR-2][i]=$i}
        ARGIND==3 && $3>0{ # ensure the actual count is nonzero
        orbit=$2; canon=orbit2canon[orbit]; edges=canonEdges[canon]; if(edges<minEdges) next;
        if(onlyBestOrbit) orbit=0;
        Kc[$1][orbit]+=$3; # increment the cluster count for appropriate orbit
        T[$1]+=$3; # keep total count of *all* orbits
        for(j=4;j<=NF;j++){ # saving the neighbors of those cliques that have high edge density for BFS
            ++neighbors[$1][orbit][$j];
        }
        }
        END{
        for(u in Kc) for(orbit in Kc[u]) {
            ORS=" "
            if('$WEIGHTED') print Kc[u][orbit]^2/T[u], u, orbit
            else            print Kc[u][orbit], u, orbit # print the near-clique count and the node
            for (v in neighbors[u][orbit]){
            print v
            }
            ORS="\n"; print "";
        }
        }' canon_maps/canon_list$k.txt canon_maps/orbit_map$k.txt $TMPDIR/blant$k.out |
    sort -gr | # > $TMPDIR/cliqs$k.sorted  # sorted near-clique-counts of all the nodes, largest-to-smallest
    hawk 'BEGIN{Srand();OFS="\t"; ID=0;}
        ARGIND==1{++degree[$1];++degree[$2];edge[$1][$2]=edge[$2][$1]=1} # get the edge list
        ARGIND==2 && !($2 in count){
        orbit=$3; count[$2]=$1; node[FNR]=$2; line[$2]=FNR;
        for(i=4; i<=NF; i++) neighbors[$2][$i] = neighbors[$i][$2]=1;
        }
        function EdgeCount(v,       edgeHits,u) {
        edgeHits=0;
        for(u in S){if(edge[u][v]) ++edgeHits;}
        return edgeHits;
        }
        function highRelCliqueCount(u, v){ # Heuristic
        if(!(u in count) || count[u]==0) return 1;
        if(!(v in count) || count[v]==0) return 0;
        if (v in count){
            return count[v]/count[u]>=0.5;
        } else {
            return 1/count[u]>=0.5;
        }
        }

        function expand(u, origin,    v,oldOrder){
        if(u in neighbors) {
            oldOrder=PROCINFO["sorted_in"];
            PROCINFO["sorted_in"]="randsort";
            for (v in neighbors[u]){
            if(!(v in visited) && (!(v in line) || line[v] > line[origin]) && highRelCliqueCount(u, v)){
                QueueAdd("Q", v);
                visited[v]=1;
            }
            }
            PROCINFO["sorted_in"]=oldOrder;
        }
        return;
        }
        END{n=length(degree); # number of nodes in the input network
        cluster[0]=1; delete cluster[0]; # cluster is now explicitly an array, but with zero elements
        numCliques=0;
        QueueAlloc("Q");
        for(start=1; start<=int(FNR*'$edgeDensity'); start++) { # look for a cluster starting on line "start". 
            if(QueueLength("Q")>0){QueueDelloc("Q");QueueAlloc("Q");} #Ensure the queue is empty
            origin=node[start];
            delete S; # this will contain the nodes in the current cluster
            delete visited;
            misses=0; # how many nodes have been skipped because they did not work?
            QueueAdd("Q", origin);
            visited[origin] = 1;
            edgeCount = 0;
            while(QueueLength("Q")>0){
            u = QueueNext("Q");
            newEdgeHits = EdgeCount(u);
            edgeCount += newEdgeHits;
            S[u]=1;
            Slen = length(S);
            if (Slen>1){
                maxEdges = choose(Slen,2);
                if(edgeCount/maxEdges < '$edgeDensity') {
                if(length(S)==Slen){
                    delete S[u]; # no node was removed, so remove this one
                    edgeCount -= newEdgeHits;
                }
                if(++misses > MAX(Slen, n/100)) break; # 1% of number of nodes is a heuristic...
                # keep going until count decreases significantly; duplicate cluster removed in the next awk
                #visited[u]=0;
                }else{
                    expand(u, orign);
                }
            } else {
                expand(u, origin);
            }
            }
            if(length(S)>'$k') {
            maxEdges=choose(length(S),2);
            ++numCliques; printf "%d %d", length(S), edgeCount
            for(u in S) {cluster[numCliques][u]=1; printf " %s", u}
            print ""
            if('$ONLY_ONE') exit;
            }
        }
        }' "$net" - | # $TMPDIR/cliqs$k.sorted | # dash is the output of the above pipe (sorted near-clique-counts)
    sort -nr | # sort the above output by number of nodes in the near-clique
    hawk 'BEGIN{ numCliques=0 } # post-process to remove duplicates
        {
        delete S;
        numNodes=$1
        edgeHits=$2;
        for(i=3;i<=NF;i++) ++S[$i]
        if(length(S)!=numNodes) next;#ASSERT(length(S)==numNodes,"mismatch in numNodes and length(S)");
        add=1;
        for(i=1;i<=numCliques;i++) {
            same=0;
            for(u in S) if(u in cluster[i])++same;
            if(same > length(cluster[i])*'$t'){add=0; break;}
        }
        if(numCliques==0 || add==1) {
                maxEdges=choose(length(S),2);
                ++numCliques; 
                printf "%d %d '$k'",length(S),edgeHits
                for(u in S) {cluster[numCliques][u]=1; printf " %s", u}
                print ""
        }
        }' | sort -k 1nr -k 4n > $TMPDIR/subfinal$k$edgeDensity.out &  # sort by number of nodes and then by the first node in the list
    done
    for k in "${Ks[@]}"; do
        wait; (( BLANT_EXIT_CODE += $? ))
    done
    sort -k 1nr -k 4n $TMPDIR/subfinal?$edgeDensity.out |
    hawk 'BEGIN{ numCliques=0 } # post-process to remove duplicates
        {
        delete S; 
        numNodes=$1
        edgeHits=$2;
        k=$3;
        for(i=4;i<=NF;i++) ++S[$i]
        if(length(S)!=numNodes) next; #ASSERT(length(S)==numNodes,"mismatch in numNodes and length(S)");
        add=1;
        for(i=1;i<=numCliques;i++) {
                same=0;
                for(u in S) if(u in cluster[i])++same;
                if(same > length(cluster[i])*'$t'){ add=0; break;}
        }
        if(numCliques==0 || add==1) {
                maxEdges=choose(length(S),2);
                ++numCliques; 
                printf "%d nodes, %d of %d edges from k %d (%g%%):",
            length(S), edgeHits, maxEdges, k, 100*edgeHits/maxEdges
                for(u in S) {cluster[numCliques][u]=1; printf " %s", u}
                print ""
        }
        }' | sort -k 1nr -k 11n > $TMPDIR/final$edgeDensity.out # sort by number of nodes and then by the first node in the list
exit $BLANT_EXIT_CODE