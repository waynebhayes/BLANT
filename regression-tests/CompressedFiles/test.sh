#!/bin/bash
die() { echo "$@" >&2; exit 1
}

echo 'Testing measurement CompressedFiles'

TEST_DIR=`pwd`/regression-tests/CompressedFiles
[ -d "$TEST_DIR" ] || die "should be run from top-level directory of the BLANT repo"

declare -a arr=("MMusculus.el.gz" "HSapiens.el.bz2" "AThaliana.el.xz")

exitCode=0

for i in "${arr[@]}"
do
   echo "Testing compressed edglist $i"
   lines=`./blant -k 3 $TEST_DIR/$i | wc -l | awk '{print $1}'`;
   if [ $lines -lt 100 ]
   then
       echo "Test 1: Failed to load network $i"
       exitCode=1
   fi
done

echo 'Done testing compression'
exit $exitCode

