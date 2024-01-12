#!/bin/bash
# Test to ensure the frequency of graphlets is within 6 sigma of the theoretical Poisson frequency based on a 3e9
# sample pre-generated "correct" frequency.
# The -t option tests parallelism, attemting to run multiple threads simultaneously.
# For now, NBE is the only sampling method that seems to get the bias bang on.

N=3000000
S=MCMC
for k in 3 4 5 6 7 8; do
    if [ -f canon_maps/canon_map$k.bin ]; then
	echo "Checking $S's raw frequency of $k-graphlets is Poisson using $N samples from networks/syeast.el"
	./scripts/test-raw-counts.sh ./blant $S $N $k regression-tests/0-sanity/syeast.$S.raw.3e9.k$k.txt |
	    awk '$1>6{printf "6-sigma violation: %s\n", $0; exit 1}' || exit 1
    fi
done
