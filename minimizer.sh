#!/bin/sh
CANON_DIR=/var/preserve/Graphette/canon_maps
hawk '
    index(FILENAME,"canon_list")&&FNR>1{# D[k][c]=decimal representation of k-c canonical; connected=Boolean
	F=FILENAME;gsub("[^0-9]"," ",F);k=1*F;D[k][FNR-2]=$1;connected[k][FNR-2]=$2;next}
    index(FILENAME,"subcanon_map"){
	F=FILENAME;gsub("[^0-9]"," ",F);k=1*F;
	    # M1 = minimizer one-step down in k; Mk[k1][k2]=minimizer from k1 down to k2, where k2 can be more than 1 step down
	    M1[k][FNR-1] = $2
	    #Mk[k][FNR-1] = ( k==4 ? $2 : M[k-1][$2] ) # k=4 is base case, otherwise sub-minimizer
	    for(i=2;i<=NF;i++){
		subcanon[k][FNR-1][i-2]=$i
		if(connected[k-1][$i]){
		    if(!connected[k-1][M1[k][FNR-1]] || $i < M1[k][FNR-1])
			M1[k][FNR-1] = $i
		}
	    }
	next
    }
    connected[$1][$2]{m1=M1[$1][$2];m2=M1[$1-1][m1];m3=M1[$1-2][m2];m4=M1[$1-3][m3]; print $1,$2,m1,m2,m3,m4}
    #END{ for(k=4;k<=8;k++) { for(c in subcanon[k]) { for(i in subcanon[k][c]) print k,c,i,subcanon[k][c][i] } } }
    ' $CANON_DIR/canon_list?.txt $CANON_DIR/subcanon_map?-?.txt "$@"
