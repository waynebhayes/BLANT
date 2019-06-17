#!/bin/bash

echo "Testing WindowRep Freq Mode"

TEST_DIR=`pwd`/regression-tests/windowRepTest
if ! [ -d $TEST_DIR ]; then
    echo "Not at top-level directory of BLANT repo" >&2
    exit 1
fi

#Can change the parameter later
w=20
k=7
n=5
sampleName="NBE"
declare -a methodName=("MIN" "MAX" "DMIN" "DMAX" "LFMIN" "LFMAX")
declare -a fnames=("SCerevisiae.el" "syeast.el" "IIDhuman.el")

numUniqRep=0
for m in "${methodName[@]}"
do
    for f in "${fnames[@]}"
    do
        cmd="./blant -k$k -w$w -p$m -s$sampleName -mf -n$n $TEST_DIR/$f"
        numUniqRep=`$cmd | sort -n  | cut -d" " -f1 | uniq | awk '{if($0 != 0) print $0}' | wc -l`
        if [ $numUniqRep -gt $n ]; then
            echo "Error($exitcode) Number of uniq windowRepInt is more than window samples. Impossible." >&2
            echo "Cmd: $cmd" >&2
            exit 1
        fi
    done
done
echo "Done Testing windowRepWindow Freq Mode"
exit 0

