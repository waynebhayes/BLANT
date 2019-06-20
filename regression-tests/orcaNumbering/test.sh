#!/bin/bash

MYTMPDIR=$(mktemp -d)
trap "rm -rf $MYTMPDIR" EXIT

for k in 3 4 5; do
    echo "Testing ODV output matches ORCA for k=${k}"
    $BLANT_DIR/blant -mo -k $k -s NBE -n 100000 $BLANT_DIR/regression-tests/orcaNumbering/$k.el | sort -n > $MYTMPDIR/$k.odv
done

python3 regression-tests/orcaNumbering/test.py $MYTMPDIR
if [[ $? -eq 0 ]]; then
    exit 0
else
    exit 1
fi
