#!/bin/sh
# Near the start of a run, we need smaller samples (say 10,000) and shorter stagnation (eg, 1000).
# As things improve, we can up the samples to 100,000 and stagnation to eg 10,000.
k=5
S=100000
STAG=1000
while true; do
    ./blant -k $k -mi -s $S syeast.random > blant.ryeast.k$k
    ./blant -k $k -mi -s $S networks/syeast.el > blant.syeast.k$k
    make synthetic
    md5sum syeast.random
    ./synthetic -s $STAG -k $k networks/syeast.el syeast.random blant.syeast.k$k blant.ryeast.k$k > x$$
    mv x$$ syeast.random
done
