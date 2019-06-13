#!/bin/sh
USAGE="$0 [-P t test.el ] k source.el; -P means make predictions with threshold t against test.el"

PREDICT=0
thresh=0
case "$1" in
-P) PREDICT=1; thresh="$2"; testNet="$3"; shift 3;;
esac

k=$1
sourceNet=$2
n=100000
while true; do
    echo --- $n ---
    blant -t 1 -s EBE -mj -n $n -k $k $sourceNet | hawk '
    ARGIND==1{edge[$1][$2]=edge[$2][$1]=1;D[$1]++;D[$2]++}
    ARGIND==2&&$2==1{numEdges[FNR-2]=$3}
    ARGIND==3{EHD[$1][$2]=EHD[$2][$1]=$3}
    ENDFILE{if(ARGIND==3){
	#EHDmap[6][8][2][2]=1 # g1 g2 o1 o2
	#EHDmap[6][8][3][2]=1
	six[6]=1;
	for(i=1;i<'$k'^2/2;i++)for(g1 in numEdges)for(g2 in numEdges)if(EHD[g1][g2]==i&&numEdges[g2]==numEdges[g1]+i)
	    EHDmap[g1][g2][2][2]=i # g1 is i edges away from g2 and we intialize both the 2nd columns with this fact.
	    # We will update the other columns below since we do not yet know how many there are.
    }}
    ARGIND==4{
	g=$1
	for(o=2;o<=NF;o++){
	    #The next two lines update the number of columns (ie., orbits) for each graphlet, up to NF
	    if (g  in EHDmap)for(g2 in EHDmap[g])EHDmap[g][g2][o][2]=EHDmap[g][g2][2][2];
	    for(g1 in EHDmap) if(g  in EHDmap[g1])for(o1 in EHDmap[g1][g])EHDmap[g1][g][o1][o]=EHDmap[g1][g][2][2];
	    n=split($o,a,":")
	    # Update the frequency that the all nodes-pairs in orbit o co-uccur in orbit o of graphlet g.
	    for(i=1;i<n;i++)for(j=i+1;j<=n;j++){
		u=MIN(a[i],a[j]); v=MAX(a[i],a[j])
		K[g][u][v][o]++;
	    }
	}
    }
    END{
	for(g1 in EHDmap)for(g2 in EHDmap[g1])for(o1 in EHDmap[g1][g2])for(o2 in EHDmap[g1][g2][o1]){
	    id=g1"("o1")"g2"("o2")";IDs[id]=1
	    # Compute the correlation that (u,v) occur as an edge in g2 and not an edge in g1.
	    # In other words, the resulting correlation of this id will tell us the following:
	    # if you are a pair of *unconnected* nodes (u,v) sharing orbit o1 in g1,
	    # how likely is that you should actually *be* attached because adding the edge (u,v)
	    # makes you look like the much more frequent (u,v) edge at orbit o2 in g2?
	    for(u in K[g2])for(v in K[g2][u])if( edge[u][v])PearsonAddSample(id,K[g2][u][v][o2],1)
	    for(u in K[g1])for(v in K[g1][u])if(!edge[u][v])PearsonAddSample(id,K[g1][u][v][o1],0)
	}
	for(id in IDs){PearsonCompute(id);if(_Pearson_rho[id])printf "%s %s\n",id,PearsonPrint(id)}
    }' $sourceNet canon_maps/canon_list$k.txt canon_maps/EdgeHammingDistance$k.txt - | sort
    n=`expr $n '*' 10`
done
