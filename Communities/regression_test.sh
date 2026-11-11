#!/bin/bash

# Currently only supporting one edgelist at a time
if [ $# -eq 0 ]; then
    echo "Please specify file to run on as an argument."
    echo "Provide an edge list to run tests on OR a test file."
    echo "The test file should be formatted as: NETWORK OBJECTIVE_FUNCTION EXPECTED_SCORE ACCEPTABLE_THRESHOLD"
    echo "You may optionally specify:"
    echo "			      1) -s: Set a scoring function. Default is HayesScore (refer to measures.c to see what is implemented)"
    echo "			      2) -n: Set the number of tests to run with. Default is 10 tests."
    echo "			      3) -e: Set an expected score."
    echo "			      4) -t: Set this flag to treat it as a test file instead of an edge list"
    exit 
fi

PASS=0
FAIL=0
EXPECTED=""
NUMRUNS=10
EDGELIST=$1
SCORINGFUNC="HayesScore"
TESTFILEMODE=0
ERRS=0
shift
echo "Running regression tests..."

while getopts "n:s:e:t" flag; do
    case "${flag}" in
	n) echo "Running $OPTARG" times ; NUMRUNS=$OPTARG;;
	e) echo "Expected score = $OPTARG" ; EXPECTED="Final score = $OPTARG";;
	s) echo "Setting scoring function to $OPTARG"; SCORINGFUNC=$OPTARG;;
	t) echo "Running in test file mode"; TESTFILEMODE=1;;
    esac
done

if [[ $TESTFILEMODE == 1 ]]; then
    while read -r NETWORK OBJECTIVE_FUNC MEAN FOUR_SIGMA; do
	echo "Run test $NETWORK $OBJECTIVE_FUNC $MEAN $FOUR_SIGMA"
	OUTPUT=$(./a.out "$NETWORK" -s "$OBJECTIVE_FUNC" 2>/dev/null)
	ACTUAL=$(echo "$OUTPUT" | grep "Final score" | awk '{print $NF}')
	LOW=$(echo "$MEAN - $FOUR_SIGMA" | bc -l)
	HIGH=$(echo "$MEAN + $FOUR_SIGMA" | bc -l)
	if (( $(echo "$ACTUAL >= $LOW && $ACTUAL <= $HIGH" | bc -l) )); then
	    echo "$NETWORK $OBJECTIVE_FUNC $ACTUAL PASS"
	else
	    echo "FATAL: does not match $ACTUAL"
	    ((ERRS++))
	fi
    done < "$EDGELIST"
    exit $ERRS
fi



if [[ $EXPECTED != "" ]]; then
    echo "Set expected score to $2"
    for i in $(seq 1 $NUMRUNS); do
	OUTPUT=$(./a.out $EDGELIST -s $SCORINGFUNC 2>/dev/null)
	ACTUAL=$(echo "$OUTPUT" | grep "Final score")
	
	if [[ "$ACTUAL" =~ "$EXPECTED" ]]; then
	    echo "Run $i: PASS"
	    ((PASS++))
	else
	    echo "Run $i: FAIL ('$ACTUAL')"
	    ((FAIL++))
	fi
    done
    echo ""
    echo "Results: $PASS passed, $FAIL failed"
    exit 
else
    echo "No expected score set"
    declare -A scores
    for i in $(seq 1 $NUMRUNS); do
	OUTPUT=$(./a.out $EDGELIST -s $SCORINGFUNC 2>/dev/null)
	ACTUAL=$(echo "$OUTPUT" | grep "Final score")
	if [[ -v scores[$ACTUAL] ]]; then
	    ((scores[$ACTUAL]++))
	else
	    scores["$ACTUAL"]=1
	fi
    done
    for score in "${!scores[@]}"; do 
	 echo "  $score — ${scores[$score]}/$NUMRUNS runs"	
    done
fi
