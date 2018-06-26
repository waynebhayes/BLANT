#!/bin/sh
CANON_DIR=/var/preserve/Graphette/canon_maps
hawk '
    index(FILENAME,"canon_list")&&FNR>1{
	F=FILENAME;gsub("[^0-9]"," ",F);k=1*F;D[k][FNR-2]=$1;connected[k][FNR-2]=$2;next}
    index(FILENAME,"subcanon_map"){
	F=FILENAME;gsub("[^0-9]"," ",F);k=1*F;M[k][FNR-1]=$2;for(i=2;i<=NF;i++){subcanon[k][FNR-1][i-2]=$i;if(connected[k-1][$i]){if(!connected[k-1][M[k][FNR-1]] || $i < M[k][FNR-1]) M[k][FNR-1] = $i}}next}
    connected[$1][$2]{m1=M[$1][$2];m2=M[$1-1][m1];m3=M[$1-2][m2];m4=M[$1-3][m3]; print $1,$2,m1,m2,m3,m4}
    #END{ for(k=4;k<=8;k++) { for(c in subcanon[k]) { for(i in subcanon[k][c]) print k,c,i,subcanon[k][c][i] } } }
    ' $CANON_DIR/canon_list?.txt $CANON_DIR/subcanon_map?-?.txt "$@"
