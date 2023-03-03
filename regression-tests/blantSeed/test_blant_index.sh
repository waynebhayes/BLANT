#!/bin/bash
[ "$CI" = true ] && exit 0

TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
./blant -k8 -lDEG2 -sINDEX -mr -a1 networks/syeast0/syeast0.el >$TMPDIR/syeast0.index 2>/dev/null
./dedup.sh $TMPDIR/syeast0.index

echo Checking for index correctness
if diff $TMPDIR/syeast0.index seed_mining/examples/syeast0.index; then
    :
else
    echo "ERROR: indexes different; see $TMPDIR" >&2
    trap ""
    exit 1
fi
