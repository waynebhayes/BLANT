#!/bin/bash
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
 trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N
python3 seed_mining/odv_helpers.py networks/syeast0/syeast0.el networks/syeast05/syeast05.el seed_mining/examples/syeast0-k5.odv seed_mining/examples/syeast05-k5.odv 5 1004 >$TMPDIR/syeast0-syeast05-k5-n1004.odvseeds 2>/dev/null

if diff $TMPDIR/syeast0-syeast05-k5-n1004.odvseeds seed_mining/examples/syeast0-syeast05-k5-n1004.odvseeds; then
    echo "odv seeds different"
    exit 1
fi
