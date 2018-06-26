#!/bin/sh
CANON_DIR=/var/preserve/Graphette/canon_maps
hawk '
    index(FILENAME,"canon_list")&&FNR>1{# D[k][c]=decimal representation of k-c canonical; connected=Boolean
	F=FILENAME;gsub("[^0-9]"," ",F);k=1*F;D[k][FNR-2]=$1;connected[k][FNR-2]=$2;next}
    index(FILENAME,"subcanon_map"){
	F=FILENAME;gsub("[^0-9]"," ",F);k=1*F;
	    ord=FNR-1 # the ordinal of the canonical on this line.
	    # M1[k][ord]=minimizer one-step down in k
	    curM1=M1[k][ord]=$2 # current M1 minimizer---just initializing, actual minimum is searched below
	    # Mk[k][ord][k2] = k2-minimizer for k-sized (ordinal) canonical ord, where k2 can be more than 1 step below k.
	    Mk[k][ord][k] = ord # every k-canonical is its own k-minimizer.
	    for(i=2;i<=NF;i++){
		removed_node=i-2;
		subcanon[k][ord][removed_node]=$i
		if(connected[k-1][curM1]){ # need to separate the cases of connected or not
		    # replace if new one is both connected and smaller
		    if(connected[k-1][$i] && $i < curM1)
			curM1 = M1[k][ord] = Mk[k][ord][k-1] = $i
		}
		else {
		    if(connected[k-1][$i] || $i < curM1) # replace if new one is smaller, or connected
			curM1 = M1[k][ord] = Mk[k][ord][k-1] = $i
		}
	    }

	    for(k2=3; k2<k;k2++)
	    {
		curMk[k2] = Mk[k][ord][k2] = ( k==4 ? $2 : Mk[k-1][$2][k2] ) # k=4 is base case, otherwise sub-minimizer
		if(connected[k2][curMk[k2]]){
		    # note removed_node for k2 may go beyond , which is OK in awk but not other languages!
		    candidate = Mk[k][ord][k2];
		    if(connected[k2][candidate] && candidate < curMk[k2])
			curMk[k2] = Mk[k][ord][k2] = candidate;
		}
		curMk[k2] = Mk[k][ord][k2] = ( k==4 ? $2 : Mk[k-1][$2][k2] ) # k=4 is base case, otherwise sub-minimizer
	    }
	next
    }
    connected[$1][$2]{m1=M1[$1][$2];m2=M1[$1-1][m1];m3=M1[$1-2][m2];m4=M1[$1-3][m3]; print $1,$2,m1,m2,m3,m4}
    #END{ for(k=4;k<=8;k++) { for(c in subcanon[k]) { for(i in subcanon[k][c]) print k,c,i,subcanon[k][c][i] } } }
    ' $CANON_DIR/canon_list?.txt $CANON_DIR/subcanon_map?-?.txt "$@"
