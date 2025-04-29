#!/bin/sh
PATH=/home/sana/bin:$PATH
die() { echo "$@" >&2; exit 1
}
[ $# -eq 2 ] || die "Must be exactly 2 input files: the edge lists of the two networks to compare"

TMP=/tmp/pushkar$$
 trap "/bin/rm -rf $TMP" 0 1 2 3 15
mkdir -p $TMP/1 $TMP/2
echo $TMP

K_Stest() { /home/sana/bin/K-Stest "$@"
}
parse() { /home/sana/bin/parse.awk "$@"
}
awkcel() { /home/sana/bin/awkcel "$@"
}
#Eigenvalues 0:
#Avg. clustering coefficient: 0.055996418392625084
#Diameter: 7
#K-hop distribution: 1504 11558 272344 1160590 724276 83410 2320 6
#Degree distribution: 0 361 215 158 130 105 76 56 37 53 38 26 17 17 9 14 13 19 16 15 16 11 7 9 
pushkar() {
    #python3 main.py -e 0 "$1" > $TMP/$2/out
    ./main -e 0 "$1" > $TMP/$2/out
    fgrep : $TMP/$2/out | egrep -v 'Eigen|NOT' > $TMP/$2/header.txt
    sed -n '/nodeName/,/^node1 node2/p' $TMP/$2/out | sed -e 's/^[# ]*//' -e 's/ /	/g' | fgrep -v edge_betweenness > $TMP/$2/nodeVals.tsv
    sed -n '/node1 node2/,/^###/p' $TMP/$2/out | sed -e 's/^[# ]*//' -e 's/ /	/g' | fgrep -v End-of-output > $TMP/$2/edgeVals.tsv
    awk '/Degree distribution:/{for(i=3;i<=NF;i++)for(c=0;c<$i;c++)print i-3}' $TMP/$2/out > $TMP/$2/deg-dist
    awk '/K-hop distribution:/{for(i=3;i<=NF;i++)for(c=0;c<$i/100;c++)print i-3}' $TMP/$2/out > $TMP/$2/k-hop
}

pushkar "$1" 1
pushkar "$2" 2

for nodeVal in clusCoff eccentricity node_betweenness; do
    for g in 1 2; do
	awkcel "{print $nodeVal}" $TMP/$g/nodeVals.tsv > $TMP/$g/$nodeVal
    done
done
for edgeVal in edge_betweenness; do
    for g in 1 2; do
	awkcel "{print $edgeVal}" $TMP/$g/edgeVals.tsv > $TMP/$g/$edgeVal
    done
done

########################################################### Global
#Eigenvalues 0:
#Nodes: 1004 Edges: 8323
#Connected components: 1
#CC byNode: 1004
#CC byEdge: 8323
#Transitivity (3* triangles/all-triads) [NOT COMPUTED]: 0
#Avg. clustering coefficient: 0.6473689095751828
#Diameter: 15
#K-hop distribution: 16646 51938 97622 155390 183390 181682 144932 96398 49304 20390 6868 1876 450 108 18
#Degree distribution: 0 104 73 51 56 49 48 41 22 32 24 26 20 28 31 27 17 19 14 9 17 20 24 20 25 21 17 12 11 8 3 6 1 5 4 5 3 5 5 5 5 1 2 3 4 3 2 1 4 2 2 3 2 3 2 0 3 0 1 0 3 1 1 1 1 2 2 1 2 1 2 2 2 1 6 4 2 1 0 1 1 2 1 1 3 1 1 1 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1
#nodeName clusCoff eccentricity node_betweenness
paste $TMP/?/header.txt |
	awk '/Nodes:/{print $1,$3,$2,$4,$6,$8}
		/^Connected/{print $1,$2,$3,$6}
		/^Avg. clust/{printf "Avg. ClustCoff %g %g\n", $4,$8}
		/^Diam/{print $1,$2,$4}'
for val in 'CC byNode' 'CC byEdge' 'K-hop' Degree; do
    for i in 1 2; do
	grep "^$val" $TMP/$i/header.txt
    done
done
for val in clusCoff eccentricity node_betweenness edge_betweenness deg-dist k-hop; do
    echo -n "$val "
	N1=`wc -l < $TMP/1/$val`
	N2=`wc -l < $TMP/2/$val`
	m=`parse "MIN($N1,$N2)"`
	M=`parse "MAX($N1,$N2)"`
    K_Stest $TMP/?/$val
done
