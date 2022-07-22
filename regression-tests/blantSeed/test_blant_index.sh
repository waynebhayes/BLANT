#!/bin/bash
TEST_DIR=$(pwd)/regression-tests/blantSeed
./blant -k8 -lDEG2 -sINDEX -mr -a1 networks/syeast0/syeast0.el >$TEST_DIR/syeast0.index 2>/dev/null
./dedup.sh $TEST_DIR/syeast0.index

DIFF=$(diff $TEST_DIR/syeast0.index examples/syeast0.index)

if [ "$DIFF" == "" ]; then
    echo "indexes identical"
    exit 0
else
    echo "indexes different"
    echo "diff is $DIFF"
    exit 1
fi
