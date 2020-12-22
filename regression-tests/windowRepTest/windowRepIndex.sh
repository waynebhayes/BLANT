#!/bin/bash

echo "Testing WindowRep Index Mode"

TEST_DIR=`pwd`/regression-tests/windowRepTest
if ! [ -d $TEST_DIR ]; then
    echo "Not at top-level directory of BLANT repo" >&2
    exit 1
fi

remove_temp_file () {
    rm -rf $TEST_DIR/.regression_test_out.txt
    rm -rf $TEST_DIR/.regression_test_error.txt
}

# Can adjust the parameter later
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

for s in "${sampleName[@]}"
do
    for m in "${methodName[@]}"
    do
        for f in "${fnames[@]}"
        do
            if ! [ -f networks/$f ]; then
                echo "Cannot find $f in the networks folder." >&2
                exit 1
            fi

            cmd="./blant -k$k -w$w -p$m -M0 -s$s -mi -n$n networks/$f"
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

            cat $TEST_DIR/.regression_test_out.txt | awk -v w="$w" -v k="$k" -v n="$n" '
                    BEGIN{numWindowRep=0; numWindowRepCounter=0; numWindow=0;}
                    {
                        if(NF==w ) {
                            numWindow+=1;
                            if(numWindowRep < numWindowRepCounter) {
                                print "Error: Output numWindowRep", numWindowRepCounter, "is larger than found", numWindowRepCounter;
                                exit 1;
                            }
                            numWindowRepCounter = 0;
                        } else if (NF == 2) {
                            numWindowRep=$2;
                        } else if (NF == k) {
                            numWindowRepCounter += 1;
                        }
                    }
                    END { if (numWindow != n) exit 1;}
                '
            exitcode=$?
            if [ $exitcode -ne 0 ]; then
                echo "Error($exitcode) from cmd: $cmd" >&2
                remove_temp_file
                exit $exitcode
            fi
        done
    done
done
remove_temp_file
echo "Done Testing windowRepWindow Index Mode"
exit 0

