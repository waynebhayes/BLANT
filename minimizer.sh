#!/bin/sh
CANON_DIR=/var/preserve/Graphette/canon_maps
hawk 'index(FILENAME,"subcanon_map")>0&&NF>0{F=FILENAME;gsub("[^0-9]"," ",F);k=1*F;D[k][FNR-1]=$1;M[k][FNR-1]=$2;for(i=2;i<=NF;i++)subcanon[k][FNR-1][i-2]=$i;next}
    {m1=M[$1][$2];m2=M[$1-1][m1];m3=M[$1-2][m2]; print m1, m2, m3}
    #END{ for(k=4;k<=8;k++) { for(c in subcanon[k]) { for(i in subcanon[k][c]) print k,c,i,subcanon[k][c][i] } } }
    ' $CANON_DIR/subcanon_map*.txt "$@"
