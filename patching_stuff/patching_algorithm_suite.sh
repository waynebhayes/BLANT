#!/bin/bash

species1=$1
species2=$2

echo running for $species1 and $species2

index_dirpath="./indexes/"
s1_filename="${index_dirpath}${species1}_k8_l2_ninf_oldBlant.txt"
s2_filename="${index_dirpath}${species2}_k8_l2_ninf_oldBlant.txt"

for NUM_MATCHING in 3 4 5 6 7
do
    for PROX in 1 4
    do
        python3 patching_algorithm.py mouse rat $s1_filename $s2_filename $NUM_MATCHING $PROX $PROX False
    done
done

