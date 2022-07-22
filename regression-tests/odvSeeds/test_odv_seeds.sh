#!/bin/bash
TEST_DIR=$(pwd)/regression-tests/odvSeeds
python3 seed_mining/odv_helpers.py networks/syeast0/syeast0.el networks/syeast05/syeast05.el examples/syeast0-k5.odv examples/syeast05-k5.odv 5 1004 >$TEST_DIR/syeast0-syeast05-k5-n1004.odvseeds 2>/dev/null

DIFF=$(diff $TEST_DIR/syeast0-syeast05-k5-n1004.odvseeds examples/syeast0-syeast05-k5-n1004.odvseeds)

if [ "$DIFF" == "" ]; then
    echo "odv seeds identical"
    exit 0
else
    echo "odv seeds different"
    echo "diff is $DIFF"
    exit 1
fi
