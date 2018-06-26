#!/bin/sh
CANON_DIR=/var/preserve/Graphette/canon_maps
awk 'index(FILENAME,"subcanon_map")>0&&NF>0{F=FILENAME;gsub("[^0-9]"," ",F);k=1*F; for(i=0;i<NF;i++)subcanon[k][FNR-1][i]=$(i+1)}
    END{
	for(k=4;k<=8;k++)
	{
	    for(c in subcanon[k])
	    {
		for(i in subcanon[k][c])
		    print k,c,i,subcanon[k][c][i]
	    }
	}
    }
    ' $CANON_DIR/subcanon_map*.txt
