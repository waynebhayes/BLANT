#!/bin/bash

echo "Testing recursive Dijkstra"

exitCode=0

echo "=====Generating 10 alignments====="
module load python/3.6.8
./Dijkstracmd

echo "=====Comparing .log with past log file====="
python3 checkDiff.py oldlog.log seed7/IIDmouse_IIDhuman_7.log 
if [ $? != 0 ];
then 
    echo "Failed to pass!"
    exitCode=1
fi

rm -rf seed7/

echo "Done testing recursive Dijkstra"

exit $exitCode
