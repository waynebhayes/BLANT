#!/bin/bash

# Written by Ethan Chennault
# Intended to tests if the multithreading produces consistent result given the same seed.
# This is an issue because each thread must have a seperate seed in order to independently sample from each other
# However, despite each having a seperate seed, if the base seed specified by -r is the same
# Output must still be exactly the same.

# Must work for all output modes and sampling methods and all K values


SEED=1
N=10000
OUT_FILE="/dev/null"
GRAPH_FILE="./../networks/syeast.el"

OUTPUT_METHODS=(f o g i j)
SAMPLING_METHODS=(NBE EBE MCMC)
AVAIL_THREADS=$(getconf _NPROCESSORS_ONLN)
THREAD_COUNTS=(1 $AVAIL_THREADS)
K_VALUES=$(seq 4 7)

TEST_COUNTER=1
EXIT_CODE=0

# do -q to do quiet mode, where only show output results if they differ
# do -v to do verbose mode, showing output results regardless
QUIET=true
while getopts ":q:v" opt; do
    case $opt in
        q)
            QUIET=true
            ;;
        v)
            QUIET=false
            ;;
        *)
            ;;
    esac
done

> "$OUT_FILE"
echo "Running Thread Seeding Test. This test runs BLANT commands twice, sorts the output, and ensures they are the exact same."

for OUTPUT_METHOD in "${OUTPUT_METHODS[@]}"; do
  echo "Testing for output method $OUTPUT_METHOD (utilizing up to $AVAIL_THREADS threads) and comparing output..."
  for SAMPLING_METHOD in "${SAMPLING_METHODS[@]}"; do
    for K in $K_VALUES; do
      for NUM_THREADS in "${THREAD_COUNTS[@]}"; do
        TMP1=$(mktemp)
        TMP2=$(mktemp)

        # COMPARISON="compare.txt"
        
        DIFFERENT=0

        CMD="./blant -r $SEED -m $OUTPUT_METHOD -s $SAMPLING_METHOD -t $NUM_THREADS -n $N -k $K $GRAPH_FILE"

        # we must sort the outputs because when multithreading for indexing modes,
        # the program will produce all the same results, but in a different order due to race conditions
        $CMD 2>/dev/null | sort > "$TMP1"
        $CMD 2>/dev/null | sort > "$TMP2"


        {
            echo "================== TEST #$TEST_COUNTER =================="
            echo "Running: $CMD, with seed: $SEED, twice."
            if diff -q "$TMP1" "$TMP2" > /dev/null; then
                echo "OUTPUT MATCHES"
            else
                echo "!!! OUTPUT DOES NOT MATCH !!!"
                ((DIFFERENT = 1))
            fi
        } >> "$OUT_FILE"


        if [ "$DIFFERENT" = 1 ]; then
            echo "Encountered discrepancies in results when running $CMD twice and sorting results. Please check it."
            ((EXIT_CODE = 1))
        fi

        if [ "$QUIET" = false ] || [ "$DIFFERENT" = 1 ]; then
            {
                echo "--- SHOWING RUN OUTPUT ---"
                echo "--- RUN 1 OUTPUT ---"
                cat "$TMP1"
                echo
                echo "--- RUN 2 OUTPUT ---"
                cat "$TMP2"
                echo

                # (cat "$TMP1") >> "$COMPARISON"
                # echo "-----------" >> "$COMPARISON"
                # (cat "$TMP2") >> "$COMPARISON"
                # exit 1
            } >> "$OUT_FILE"
        fi

        # > "$COMPARISON"

        rm "$TMP1" "$TMP2"

        ((TEST_COUNTER++))
      done
    done
  done
done


if [ "$EXIT_CODE" = 0 ]; then
    echo "No discrepancies detected between any of the tests. Multithread seeding works as appropriate."
else
    echo "Discrepancies detected between some of the tests. Multithread seeding DOES NOT work correctly; consistent output should be produced with seeding."
fi
exit $EXIT_CODE