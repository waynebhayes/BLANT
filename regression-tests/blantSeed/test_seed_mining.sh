#!/bin/bash
TEST_DIR=$(pwd)/regression-tests/blantSeed
cd seed_mining
./full_algorithm_helpers.py syeast0 syeast05 ../examples/syeast0.index ../networks/syeast0/syeast0.el ../examples/syeast05.index ../networks/syeast05/syeast05.el >$TEST_DIR/syeast0-syeast05-results.txt
cd ..

DIFF=$(diff $TEST_DIR/syeast0-syeast05-results.txt examples/syeast0-syeast05-results.txt)

if [ "$DIFF" == "" ]; then
    echo "seed + extract results identical"
    exit 0
else
    echo "seed + extract results different"
    echo "diff is $DIFF"
    exit 1
fi

