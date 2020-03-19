#!/bin/sh
NUM_FAIL=0
for dir in regression-tests/*; do
    [ -d "$dir" ] || continue;
    echo --- in directory $dir ---
    for r in $dir/*.sh; do
	echo --- running test $r ---
	if "$r"; then
	    :
	else
	    NUM_FAIL=`expr $NUM_FAIL + 1`
	fi
    done
done
exit $NUM_FAIL
