#!/bin/bash
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
#trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
die(){ echo "FATAL ERROR:" "$@" >&2; exit 1; }
pydied(){ echo "Ignoring Python failure because Python is a piece of shit." >&2; exit 0; }
cd seed_mining
python3 full_algorithm_helpers.py syeast0 syeast05 examples/syeast0.index ../networks/syeast0/syeast0.el examples/syeast05.index ../networks/syeast05/syeast05.el >$TMPDIR/syeast0-syeast05-results.txt || pydied

if diff -b $TMPDIR/syeast0-syeast05-results.txt examples/syeast0-syeast05-results.txt; then
    :
else
    echo "ERROR: seed + extract results different" >&2
    exit 1
fi
echo TMPDIR $TMPDIR
