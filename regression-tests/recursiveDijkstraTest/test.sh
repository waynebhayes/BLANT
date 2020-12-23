#!/bin/bash
cd `dirname $0`
echo "Testing recursive Dijkstra"

exitCode=0

rm -rf /tmp/seed7/
mkdir -p /tmp/seed7; ln -sf /tmp/seed7 .

echo "=====Generating 10 alignments====="
if module avail 2>/dev/null; then
    module unload python
    module load python/3.6.8
fi
./Dijkstracmd
rm -f *.pickle

echo "=====Comparing .log with past log file====="
for i in oldlog.log seed7/*.log; do echo `sed 's/time:[0-9.:]*//' $i | sort | md5sum` $i; done | 
    awk 'BEGIN{FAIL=1}{print NR, $0; md5[NR]=$1}END{exit( md5[1]!=md5[2] ? FAIL : !FAIL)}'
exit $?
