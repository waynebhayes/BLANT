#!/bin/sh

# compile the file src/blant-predict-release.c using libwayne libraries:
if [ $# -ne 1 ]; then echo "must provide number of predictions as first argument" >&2; exit 1; fi

TOP="$1" # number of predictions to take
LIBWAYNE_HOME=`cd libwayne && /bin/pwd`
PATH="$PATH:$LIBWAYNE_HOME/bin:$LIBWAYNE_HOME/../scripts"
export LIBWAYNE_HOME PATH
# compile the executable
[ -f libwayne/libwayne.a ] || (cd libwayne && make all)
wgcc -O3 -Isrc -o blant-predict-release src/blant-predict-release.c

# create a random-order list of HI-union edges
TMP=/tmp/HI-union
randomizeLines < networks/HI-union.el > $TMP.el

# separate into 90% train, 10% test
FOLD=`wc -l < $TMP.el` # number of edges in network
FOLD=`expr $FOLD / 10` # take 10% for testing
head -$FOLD $TMP.el > $TMP.test
FOLD1=`expr $FOLD + 1` # training size
tail -n +$FOLD1 $TMP.el > $TMP.train

# Now run predictions
./blant-predict-release -k4 -mp6,0,0 $TMP.train | # estimate L3-paths on 4-node graphlets as in Kovacs et al.
    sort -gr | # sort the output node-pairs by L3 count, highest-to-lowest
    cut -f2 | # extract the node-pairs
    hawk '{OFS="\t"; print MAX($1,$2), MIN($1,$2)}' | # print the tab-separated node-pairs in correct direction
    fgrep -vf $TMP.train | # remove the edges that are in the training set
    head -$TOP | # take the top 1,000 predicted edges
    fgrep -xcf - $TMP.test # count the correct
