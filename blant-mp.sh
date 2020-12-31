#!/bin/sh
USAGE="USAGE: $0 [-async] blant-exec {blant args without specifying mode, which must be -mp}"

ASYNC=false
case "$1" in
-a*) ASYNC=true;shift;;
esac

BLANT=$1; shift
k=`echo "$@" | sed 's/.*-k//' | awk '{print $1}'`

#output looks like:
#P  u:v    e  o:p   q:r   x:y
#P 680:265 0 37:37 38:39 517:476
#P 680:265 0 37:37 38:39 517:25
#P 265:517 1 37:38 37:39 680:476

$BLANT -mp "$@" | hawk '{
	sub("^","'$k':",$4);
	G[$2]=$3;
	#PG4[$2][$4][$5][$6]=1 # worst in terms of memory usage
	#PG3[$2][$4][$6]=1 # PG=PredictGraph, maybe q:r is irrelevant, just unique edges?
	++PG[$2][$4]  # Simplest: literally just count (u:v,o:p) orbit-pair frequencies
    }
    END {
	for(uv in PG) {
	    printf "%s %d",uv,G[uv]
	    for(op in PG[uv]) {
		#TP4 = 0 # TP=TotalPredict of edge-orbit-pair (o,p) on node pair (u,v) in G
		#for(qr in PG4[uv][op]) TP4 += length(PG4[uv][op][qr]) # all the 1s of the xys
		#printf "\t%s %d",op,length(PG3[uv][op])
		printf "\t%s %d",op,PG[uv][op]
		#if(length(PG[uv][op])!=TP4)printf("ERROR");
	    }
	    print ""
	}
    }' &
if $ASYNC; then
    pid=`jobs -l | head -1 | awk '{print $2}'`

    while free -g | awk '/Mem:/&&$NF<2{print;++flag}/Swap:/&&$NF<2{print;++flag}flag==2{exit(flag)}'; do
	sleep 10; free -g|head -1;
    done; kill $pid # nuke the $BLANT generating the stuff for hawk
fi
wait
