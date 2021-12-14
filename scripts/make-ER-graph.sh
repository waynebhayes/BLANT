#!/bin/sh
N=`count.el "$1" | awk '{print $1}'`
cat "$1" | hawk 'BEGIN{srand()}END{while(n<NR){v=u=int('$N'*rand());while(v==u)v=int('$N'*rand()); x=MIN(u,v);y=MAX(u,v);if(!seen[x" "y]){n++;print x,y;seen[x" "y]=1}}}'
