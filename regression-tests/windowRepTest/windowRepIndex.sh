#!/bin/bash

echo "Testing WindowRep Index Mode"

TEST_DIR=`pwd`/regression-tests/windowRepTest
if ! [ -d $TEST_DIR ]; then
    echo "Not at top-level directory of BLANT repo" >&2
    exit 1
fi

# Can adjust the parameter later
w=20
k=7
n=5
sampleName="NBE"
declare -a methodName=("MIN" "MAX" "DMIN" "DMAX" "LFMIN" "LFMAX")
declare -a fnames=("SCerevisiae.el" "syeast.el" "IIDhuman.el")

for m in "${methodName[@]}"
do
    for f in "${fnames[@]}"
    do
        if ! [ -f networks/$f ]; then
            echo "Cannot find $f in the networks folder." >&2
            exit 1
        fi
        cmd="./blant -k$k -w$w -p$m -s$sampleName -mi -n$n networks/$f"
        $cmd| awk -v w="$w" -v k="$k" -v n="$n" '
                BEGIN{numWindowRep=0; numWindowRepCounter=0; numWindow=0;}
                {
                    if(NF==w ) {
                        numWindow+=1;
                        if(numWindowRep != numWindowRepCounter) exit 1;
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
            exit $exitcode
        fi
    done
done
echo "Done Testing windowRepWindow Index Mode"
exit 0

