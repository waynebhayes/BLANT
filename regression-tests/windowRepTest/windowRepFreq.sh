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
declare -a sampleName=("MCMC" "NBE")
declare -a methodName=("MIN" "MAX" "DMIN" "DMAX" "LFMIN" "LFMAX")
declare -a fnames=("SCerevisiae.el" "AThaliana.el" "CElegans.el")

numUniqRep=0
for s in "${sampleName[@]}"
do
    for m in "${methodName[@]}"
    do
        for f in "${fnames[@]}"
        do
            if ! [ -f networks/$f ]; then
                echo "Cannot find $f in the networks folder" >&2
                exit 1
            fi
            cmd="./blant -k$k -w$w -p$m -s$s -mf -n$n networks/$f"
            numUniqRep=`$cmd | awk '{if($1 != 0) print $2}' | wc -l`
            if [ $numUniqRep -gt $n ]; then
                echo "Error: Number of uniq windowRepInt is more than window samples. Impossible." >&2
                echo "Cmd: $cmd" >&2
                exit 1
            fi
            if [ $numUniqRep -eq 0 ] && [ $n -gt 0 ]; then
                echo "Error: Did not find any windowRepInt. Impossible." >&2
                echo "Cmd: $cmd" >&2
                exit 1
            fi
        done
    done
done
echo "Done Testing windowRepWindow Freq Mode"
exit 0

