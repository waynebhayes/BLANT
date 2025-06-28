#!/bin/sh
XARGS=`echo "$@" | awk '{if(!/Darwin/)print "-r"}'`
PRUNE_egrep=`ls -d Patrick .git PPI-predict HI-union Graphette-faye SYNTH* 2>/dev/null | awk '{printf "%s%s",bar,$1;bar="|"}END{print ""}'`
PRUNE_find=`echo "$PRUNE_egrep" | sed -e 's,^,-path \./,' -e 's,|, -o -path ./,g'`
find . '(' $PRUNE_find ')' -prune -o -name __pycache__ -o -name '*.pyc' | egrep -v "$PRUNE_egrep" | xargs $XARGS /bin/rm -rf &
