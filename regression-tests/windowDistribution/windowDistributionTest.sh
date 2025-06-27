#!/bin/bash
CORES=${CORES:=1}
echo "Testing Distribution Output Mode"

TEST_DIR=`pwd`/regression-tests/windowDistribution
if ! [ -d $TEST_DIR  ]; then
    echo "Not at top-level directory of BLANT repo" >&2
    exit 1
fi

remove_temp_file () {
    rm -rf $TEST_DIR/blant_distr1.txt
    rm -rf $TEST_DIR/blant_distr2.txt
}

# Can change the following parameters
n=1000
seed=0
declare -a fnames=("SCerevisiae.el" "AThaliana.el" "CElegans.el")

for k in {3..6}; do
    for nnw in "${fnames[@]}"
    do
	CMD="./blant -q -k$k -md -sMCMC -n$n networks/$nnw"
	outLines=`$CMD | tr ' ' '\n' | awk '{sum+=$1} END {print 1*sum}'`
	if [ $outLines -ne $n ]
	then
	    echo "Distribution table entries do not sum up to $n" >&2
	    echo "cmd was: $CMD" >&2
	    exit 1
	fi
	CMD="./blant -q -k$k -md -sMCMC -r$seed -n$n networks/$nnw"
	$CMD > $TEST_DIR/blant_distr1.txt
	$CMD > $TEST_DIR/blant_distr2.txt
	if ! [ `diff $TEST_DIR/blant_distr1.txt $TEST_DIR/blant_distr2.txt | wc -l` -eq 0 ]; then
	    echo "Distribution table differs for the same random seed $seed" >&2
	    echo "cmd was: $CMD" >&2
	    remove_temp_file
	    exit 1
	fi
    done
done

remove_temp_file
echo "Done Testing Distribution Output Mode"

