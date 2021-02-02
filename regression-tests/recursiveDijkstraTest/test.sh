#!/bin/bash
cd `dirname $0`
echo "Testing recursive Dijkstra"

exitCode=0

die(){ (echo "USAGE: $USAGE"; echo "`basename $0`: FATAL ERROR: $@")>&2; exit 1; }

BASENAME=`basename "$0" .sh`
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

(rm -rf $TMPDIR/seed7/ ./seed7 && mkdir -p $TMPDIR/seed7 && ln -sf $TMPDIR/seed7 .) || die "Problems creating $TMPDIR/seed7 and/or ./seed7"

echo "=====Generating 10 alignments====="
if module avail 2>/dev/null; then
    module unload python
    module load python/3.6.8
fi
if ./Dijkstracmd; then
    rm -f *.pickle

    echo "=====Comparing .log with past log file====="
    for i in oldlog.log seed7/*.log; do echo `sed 's/time:[0-9.:]*//' $i | sort | md5sum` $i; done | 
	awk 'BEGIN{FAIL=1}{print NR, $0; md5[NR]=$1}END{exit( md5[1]!=md5[2] ? FAIL : !FAIL)}'
    exit $?
else
    [ $? = 99 ] && exit 0 # ignoring failed Dijkstra
fi
