#!/bin/bash

echo "Testing WindowRep Freq Mode"

TEST_DIR=`pwd`/regression-tests/windowRepTest
if ! [ -d $TEST_DIR ]; then
    echo "Not at top-level directory of BLANT repo" >&2
    exit 1
fi

remove_temp_file () {
    rm -rf $TEST_DIR/.regression_test_out.txt
    rm -rf $TEST_DIR/.regression_test_error.txt
}

#Can change the parameter later
w=20
if [ -f canon_maps/canon_map7.bin ]; then
    k=7
else
    k=6
fi
n=5
max_depth_attempt=10
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

            cmd="./blant -k$k -w$w -p$m -M0 -s$s -mf -n$n networks/$f"
            `$cmd 1>$TEST_DIR/.regression_test_out.txt 2>$TEST_DIR/.regression_test_error.txt`

            depth_attempt=0
            while [ `grep 'Assertion \`depth++ < \| Assertion \`++numTries <' $TEST_DIR/.regression_test_error.txt | wc -l` -gt 0 ] && [ $depth_attempt -lt $max_depth_attempt ]
            do
                `$cmd 1>$TEST_DIR/.regression_test_out.txt 2>$TEST_DIR/.regression_test_error.txt`
                depth_attempt=$((depth_attempt + 1))
            done

            if [ $depth_attempt -eq $max_depth_attempt ]; then
                echo "Error:  Window Sampling Failed. Increase MAX_TRIES in BLANT and try again"
                echo "cmd:  " $cmd
                remove_temp_file
                exit 1
            fi

            numUniqRep=`cat $TEST_DIR/.regression_test_out.txt | awk '{if($1 != 0) print $2}' | wc -l`
            if [ $numUniqRep -gt $n ]; then
                echo "Error: Number of uniq windowRepInt is more than window samples. Impossible." >&2
                echo "Cmd: $cmd" >&2
                remove_temp_file
                exit 1
            fi
            if [ $numUniqRep -eq 0 ] && [ $n -gt 0 ]; then
                echo "Error: Did not find any windowRepInt. Impossible." >&2
                echo "Cmd: $cmd" >&2
                remove_temp_file
                exit 1
            fi
        done
    done
done
remove_temp_file
echo "Done Testing windowRepWindow Freq Mode"
exit 0

