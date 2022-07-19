#!/bin/bash
cd `dirname $0`

die(){ (echo "USAGE: $USAGE"; echo "`basename $0`: FATAL ERROR: $@")>&2; exit 1; }

[ -f ../../Dijkstra/libcalci.so ] || ../../libwayne/bin/wgcc -o ../../Dijkstra/libcalci.so -fPIC -shared ../../libwayne/src/graph.c ../../libwayne/src/sets.c ../../libwayne/src/tinygraph.c ../../libwayne/src/misc.c ../../libwayne/src/bintree.c ../../libwayne/src/combin.c ../../libwayne/src/queue.c ../../Dijkstra/dijkstra_helper.c || die "compile of libcalci.so failed"

BASENAME=`basename "$0" .sh`
[ $BASENAME == "$BASENAME" ] || die "something weird with filename in '$BASENAME'"

# Temporary Filename + Directory (both, you can use either, note they'll have different random stuff in the XXXXXX part)
TMPDIR=`mktemp -d /tmp/$BASENAME.XXXXXX`
#trap "/bin/rm -rf $TMPDIR; exit" 0 1 2 3 15 # call trap "" N to remove the trap for signal N

(/bin/rm -rf $TMPDIR/seed7/ ./seed7 && mkdir -p $TMPDIR/seed7 && ln -sf $TMPDIR/seed7 .) || die "Problems creating $TMPDIR/seed7 and/or ./seed7"

echo "=====Generating 10 alignments====="
if module avail 2>/dev/null; then
    module unload python
    module load python/3.6.8
fi

if ./Dijkstracmd; then
    echo "=====Comparing .log with past log file====="
    for i in oldlog.log seed7/*.log; do echo `sed 's/time:[0-9.:]*//' $i | sort | md5sum` $i; done | 
	if awk 'BEGIN{FAIL=1}{print NR, $0; md5[NR]=$1}END{exit( md5[1]!=md5[2] ? FAIL : !FAIL)}'; then
	    /bin/rm -rf $TMPDIR
	else
	    echo "logs differ--see $TMPDIR:"
	    diff -b <(sed 's/time:[0-9.:]*//' oldlog.log) <(sed 's/time:[0-9.:]*//' seed7/*.log)
	    exit 1
	fi
else
    D=$?
    [ $D = 99 ] && exit 0 # ignoring missing SIM file
    echo "Dijkstra failed" >&2
    exit $D
fi
