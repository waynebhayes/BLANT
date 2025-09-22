#!/bin/sh

# compile the file src/blant-predict-release.c using libwayne libraries:
if [ $# -ne 1 ]; then echo "must provide number of predictions as first argument" >&2; exit 1; fi

TOP="$1" # number of predictions to take

# compile the executable
[ -f libwayne/libwayne.a ] || (cd libwayne && make all)
./libwayne/bin/wgcc -O3 -Isrc -o blant-predict-release src/blant-predict-release.c

# Now run it
./blant-predict-release -k4 -mp6,0,0 networks/HI-union0.train | # estimate L3-paths on 4-node graphlets as in Kovacs et al.
    sort -gr | # sort the output node-pairs by L3 count, highest-to-lowest
    cut -f2 | # extract the node-pairs
    ./libwayne/bin/hawk '{OFS="\t"; print MIN($1,$2), MAX($1,$2)}' | # print the tab-separated node-pairs in correct direction
    fgrep -vf networks/HI-union0.train | # remove the edges that are in the training set
    head -$TOP | # take the top 1,000 predicted edges
    fgrep -xcf - networks/HI-union0.test # count the correct
