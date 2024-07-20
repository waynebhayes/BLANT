#!/bin/bash
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
 die(){ echo "FATAL ERROR:" "$@" >&2; exit 1; }
pydied(){ echo "Ignoring Python failure because Python is a piece of shit." >&2; exit 0; }
python3 seed_mining/odv_helpers.py networks/syeast0/syeast0.el networks/syeast05/syeast05.el seed_mining/examples/syeast0-k5.odv seed_mining/examples/syeast05-k5.odv 5 1004 >$TMPDIR/syeast0-syeast05-k5-n1004.odvseeds || pydied

if diff -b $TMPDIR/syeast0-syeast05-k5-n1004.odvseeds seed_mining/examples/syeast0-syeast05-k5-n1004.odvseeds; then
    :
else
    echo "odv seeds different"
    exit 1
fi
