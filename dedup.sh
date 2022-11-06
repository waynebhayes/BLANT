#!/bin/bash
deduped_path="${1}.deduped"
echo "before: `wc -l $1`"
awk '!visited[$0]++' $1 > $deduped_path
rm $1
mv $deduped_path $1
echo "after: `wc -l $1`"
